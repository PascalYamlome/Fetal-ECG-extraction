function qrsM=QRSdetectorF1(vadx,vdx,Fs,pth,RRts,pmQT)
% --------------------------------------------------------------------------------------------
% QRSdetectorF1.m: QRS detector
%
%   Input parameters:
%	 vadx : array of filtered absolute derivate values
%	 vdx  : array of filtered derivate values
%	 Fs   : sampling frequency
%	 pth  : threshold on derivative
%	 RRts : RR (seconds) typical of the animal
%	 pmQt : fraction of QT length for QT maks
%
%   Ouput parameters:
%   qrsM  : matrix (6 X number_of_detected_QRS) consisting of
%           a column of 6 time indeces for each detected QRS:
%           QRSref    (reference point, max signed derivative )
%           QRSonset  (derivative overcomes half threshold)
%           supThDer  (derivative overcomes threshold)
%           maxAbsDer (max absolute derivative)
%           infThDer  (derivative becomes lower than threshold)
%           QRSoffset (derivative decreses below half threshold)
%
%	Example:
%	 qrsM=QRSdetectorF1(vadx,vdx,fs,0.4,0.86,1);
%
% NOTE:
% Decreasing the "pth" threshold leads to sensibility increasing but specificity decreasing.
%
% 2002     Matlab translation of the "didactic" version of QRS detector developed in Mathcad.
% 2005     Inserted a coarse estimate of QRS onset and QRS offset.
% 2007     QRS fiducial point as max signed derivative
% 2008/11  Modified control row 108 (if(i < iinizio+QRSd  | i<ifine+nsd))
%
% --------------------------------------------------------------------------------------------
% Author: Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
% For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% --------------------------------------------------------------------------------------------
% This program is free software; you can redistribute it and/or modify it under the terms
% of the GNU General Public License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
%
% This program is distributed "as is" and "as available" in the hope that it will be useful,
% but WITHOUT ANY WARRANTY of any kind; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% --------------------------------------------------------------------------------------------

if (nargin < 3),
    error('At least 3 parameters are required'); return
end
if (nargin < 4),  pth=0.5; end
if (nargin < 5),  RRts=0.86;  end      % 0.86=human  % 0.25= mouse
if (nargin < 6),  pmQT=1;  end

sqRRts=sqrt(RRts);
QTlen=0.420*sqRRts;   % normal Qt length (RR=1 => QT=0.42s, humans)

%QRSd=round(0.15*RRtc);           % extension of the previous QRS detection is allowed in the interval (ifine, iinizio+QRSd)
QRSd=round((0.05+0.25*sqRRts)*Fs);    % extension of the previous QRS detection is allowed in the interval (ifine, iinizio+QRSd)
maskQT=round((0.07+pmQT*QTlen)*Fs);  % QT mask length (RR=1s => QT=0.42s, humans), QRS detection threshold decreases (linearly)
rthT=2.4;                   % from rthT*th if i=ifine to th if i=ifine+maskQt
% nsp=round(0.05*sqRRts*Fs);   % max number of sample before threshold crossing (increasing)
nsp=round(0.15*sqRRts*Fs);     % max number of sample before threshold crossing (increasing) 23-09-05
nsd=round(0.078*sqRRts*Fs);     % max number of sample after threshold crossing (decreasing )

RRtc=Fs*RRts;
nqappr=fix(1.2*length(vadx)/RRtc);   % approximate estimate of  number of QRSs
QRSref=zeros(1,nqappr); QRSonset=zeros(1,nqappr); inizio=zeros(1,nqappr);
vimaxd=zeros(1,nqappr); fine=zeros(1,nqappr); QRSoffset=zeros(1,nqappr);

isai=1.5*Fs;                      %  index of first sample used in inizialization
fsai=min(length(vadx), floor(60*RRtc));  % index of last sample used in inizialization
w2=fix(2*RRtc);           % wide windows containing at least one QRS
mD2=meanMaxSc(vadx(isai:fsai), w2, 1,1);  % compute the average of maximum derivatives on windows of 2s
% (1% of minima and 1% of maxima are discard)

meaD=mD2;
th=pth*meaD;
thbon=0.5*th;
thboff=0.5*th;
% --- choose derivative signum for fiducial QRS pointer
%[minD, maxD] =mimaxsc(vdx(1:nsai),nsc,nsc);
[minD, maxD] =mimaxsc(vdx(isai:fsai),1,1);
if(maxD > -minD*1.1), vsdx=vdx; else vsdx=-vdx; end
jq=1;
QRS=0;              % flag asserting a previous QRS threshold overcoming
maxd=0;
imaxd=1;
iinizio=-maskQT;
ifine=iinizio;
for i=1:length(vadx)    % main loop on derivative samples
    vadxi=vadx(i);
    if(QRS)             % inside QRS interval
        if(vadxi < th)    % absolute derivative becomes lower than threshold
            inizio(jq)=iinizio;  % save the index of last threshold crossing
            vimaxd(jq)=imaxd;    % save the index of max derivative
            ifine=i;
            fine(jq)=i;
            maxdc=min(maxd,th*4);               % bounding of derivative maximum
            thon=min((0.5*thbon+0.25*maxdc), th);   % adjust threshold for QRS onset
            thoff=min((0.5*thboff+0.25*maxdc), th); % adjust threshold for QRS offset
            QRSonset(jq)=max(iinizio-nsp,1)-1+min([find(vadx(max(iinizio-nsp,1):max(iinizio,1))>thon);max(iinizio-nsp,1)]);
            %            [maxs,imaxs]=max(vsdx(max(iinizio-nsp,1):min(ifine+nsd,length(vsdx))));
            [maxs,imaxs]=max(vsdx(max(iinizio,1):min(ifine,length(vsdx))));
            QRSref(jq)=max(iinizio,1)-1+imaxs;
            QRSoffset(jq)=ifine-1+max([find(vadx(ifine-1:min(ifine+nsd,length(vadx)))>thoff);1]);
            meaD=0.97*meaD + (1-0.97)* maxdc;   % updating the average of derivative maximum
            th= pth*meaD;                       % updating the threshold on derivative
            thbon=0.5*th;
            thboff=0.5*th;
            QRS=0;
            jq=jq+1;
        elseif(vadxi > maxd)
            imaxd=i;                        % save the index of max derivative
            maxd=vadxi;
        end
        
    elseif(vadxi > th)   % absolute derivative greater than threshold
        if(i < iinizio+QRSd  || i<ifine+nsd)        % extend previous QRS detection
            if(jq>1), jq=jq-1; end
            ifine=i;
            if(vadxi > maxd)
                imaxd=i;
                maxd=vadxi;
            end
            QRS=1;      % set again QRS flag, previous QRS was not ended
            % To avoid T wave the QRS detection threshold decreases (linearly) from "rthT*th" to "th" for
            % i>=ifine and i<ifine+maskQt
        elseif(vadxi > th*(rthT - (rthT-1)*(i-ifine)/maskQT))
            imaxd=i;
            iinizio=i;         % save index of the first threshold crossing
            maxd=vadxi;
            QRS=1;             % set flag, a new QRS is detected
        end
    end
end
if(jq>1)
    qrsM=[QRSref; QRSonset; inizio; vimaxd; fine; QRSoffset];
    qrsM(:,jq:end)=[];
else
    qrsM=[];
end

end %== function ================================================================
%
