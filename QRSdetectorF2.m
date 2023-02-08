function QRSref=QRSdetectorF2(vadx,vdx,Fs,pth,RRst,RRsi,lensi,fb,pmQT)
% --------------------------------------------------------------------------------------------
% QRSdetectorF2.m: QRS detector
%     QRS fiducial point as max signed derivative
%     It needs to start from a controlled signal interval
%
%   Input parameters:
%	 vadx : array of filtered absolute derivate values  
%	 vdx  : array of filtered derivate values  
%	 Fs   : sampling frequency  
%	 pth  : threshold on derivative
%	 RRst : RR (seconds) typical of the animal
%	 RRsi : RR (seconds) initial value
%    lensi: initial interval length (seconds)
%    fb   : forward (0) - backward (1)  flag
%	 pmQt : fraction of QT length for QT mask  (NOT used)
%
%   Ouput parameters:
%    QRSref    (reference point, max signed derivative )
%
%	Example:
%	 qrsM=QRSdetectorF2(vadx,vdx,fs,0.4,0.86,1);
%
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
    error ( 'At least 3 parameters are required' ); return
end
if (nargin < 4),  pth=0.5; end
if (nargin < 5),  RRst=0.86;  end  % 0.86=human  % 0.25= mouse
if (nargin < 6),  RRsi=RRst;  end
if (nargin < 7),  lensi=RRsi*10;  end
if (nargin < 8),  fb=0;  end        % 0=forward, 1=backward
if (nargin < 9),  pmQT=1;  end

RRtc=Fs*RRst;
RRci=Fs*RRsi;
lenci=lensi*Fs;
% sqRRst=sqrt(RRst);
% QTlen=0.420*sqRRst;   % normal Qt length (RR=1 => QT=0.42s, humans)

m=5;  mu=0.01; alp=mu/10;

QRSref=[];

%maskQT=round((0.07+pmQT*QTlen)*Fs);  % QT mask length (RR=1s => QT=0.42s, humans)

% % nsp=round(0.05*sqRRst*Fs);   % max number of sample before QRSrefj
% nsp=round(0.15*sqRRst*Fs);     % max number of sample before QRSrefj
% nsd=round(0.078*sqRRst*Fs);     % max number of sample after QRSrefj

wleftAmi=0.2; wrightAmi=0.2;
wleft= round(0.5*RRst*Fs);   wright= round(1.2*RRst*Fs);
wleftc= round(0.15*RRst*Fs); wrightc= round(0.15*RRst*Fs);

isai=1;                      %  index of first sample used in inizialization
fsai=min(length(vadx),round(max(4*RRci, lenci)));  % index of last sample used in inizialization

w2=fix(2*RRci);           % wide windows containing at least one QRS
mD2=meanMaxSc(vadx(isai:fsai), w2, 0.5,0.5);  % compute the average of maximum derivatives on windows of 2s
% (0.5% of minima and 0.5% of maxima are discard)
meaD=mD2;
th=pth*meaD;

% --- choose derivative signum for fiducial QRS pointer
[minD, maxD] =mimaxsc(vdx(isai:fsai),1,1);
if(fb)
    if(minD < -maxD*1.1), vsdx=-vdx; else vsdx=vdx; end  % backward
else
    if(maxD > -minD*1.1), vsdx=vdx; else vsdx=-vdx; end  % forward
end
j=1;
RRsm=RRsi;  RRcm=RRsm*Fs;     RRest=RRcm;
QRSrefj=1;

while QRSrefj<=length(vdx)    % main loop on derivative samples
    
    QRSrefjlmi=max(QRSrefj-wleft,1);   QRSrefjrma=min(QRSrefj+wright,length(vsdx));
    if(j>1)
        %---
        weightw=ones(wleft+wright+1,1);
        for k=wrightc:wright
            if(weightw(wleft+k)-0.02>wrightAmi), weightw(wleft+k+1)= weightw(wleft+k)- 0.002;
            else weightw(wleft+k+1)= weightw(wleft+k); end
        end
        for k=wleftc:wleft
            if(weightw(wleft-k+2)-0.02>wleftAmi), weightw(wleft-k+1)= weightw(wleft-k+2)- 0.002;
            else weightw(wleft-k+1)= weightw(wleft-k+2); end
        end
        weightw=weightw(1:QRSrefjrma-QRSrefjlmi+1);
        %----
        
        [maxs,imaxs]=max(weightw.*vsdx(QRSrefjlmi: QRSrefjrma));  % max of derivative
    else
        %QRSrefjrma=min(QRSrefj+round(RRcm),length(vsdx));
        QRSrefjrma=min(QRSrefj+wright,length(vsdx));
        weightw=ones(wleft+wright+1,1);
        for k=wrightc:wright
            if(weightw(wleft+k)-0.02>wrightAmi), weightw(wleft+k+1)= weightw(wleft+k)- 0.002;
            else weightw(wleft+k+1)= weightw(wleft+k); end
        end
        weightw=weightw(1:QRSrefjrma-QRSrefjlmi+1);
        [maxs,imaxs]=max(weightw.*vsdx(QRSrefjlmi: QRSrefjrma));  % max of derivative
    end
    
    if(maxs>th*weightw(imaxs))            % comparison of max of derivative with weighted threshold
        QRSrefj=QRSrefjlmi-1+imaxs;
        QRSref(j)=QRSrefj;                % QRS pointer saving
        
        if(j>1)
            RRcj=QRSref(j)-QRSref(j-1);
            RRcmp=RRcm;  RRczj=RRcj-RRcmp;
            RRcVz(j,1)=RRczj;
            
            RRcm=RRcmp+ 0.05* sign(RRczj)*min(abs(RRczj), 0.05*RRcmp);
            RRsm=RRcm/Fs;
            if(RRsm<RRst*.5 || RRsm>RRst*2.5), RRsm=RRst; RRcm=RRsm*Fs; end
%             sqRRsm=sqrt(RRsm);
%             QTlen=0.420*sqRRsm;
            %    maskQT=round((0.07+pmQT*QTlen)*Fs);
            %    nsp=round(0.15*sqRRsm*Fs);
            %    nsd=round(0.078*sqRRsm*Fs);
            wleft = round(0.3*RRsm*Fs);  wright = round(0.8*RRsm*Fs);
            %  wleftc= round(0.15*RRsm*Fs); wrightc= round(0.15*RRsm*Fs);
            wleftc= round(0.02*RRsm*Fs); wrightc= round(0.02*RRsm*Fs);
            RRest=RRcm;
            if(j==m+1)
                RRcVzj=RRcVz(2:end);
                a = lpc(RRcVzj,m); w=-a(end:-1:2)';
                RRczest=w'*RRcVzj;
                RRest=RRcmp+RRczest;             % Estimated RR value
                sz2=RRcVzj'*RRcVzj/m;
                %RRest = filter([0 -a(2:end)],1,RRcVzj);
            elseif(j>m+1)
                sz2 = sz2 + alp*sign(RRczj*RRczj-sz2)*min(RRczj*RRczj-sz2, sz2);
                cp=0.01*sz2;
                er=RRczj-RRczest;
                RRcVzj=RRcVz(j-m+1:j);
                if(er*er < 0.5*sz2)
                    w = w + (2*mu/(cp+RRcVzj'*RRcVzj))*er*RRcVzj;
                else
                    w=0.97*w;
                end
                RRczest=w'*RRcVzj;
                if(abs(RRczest)>0.03*RRcmp)
                    RRczest=sign(RRczest)*0.03*RRcmp;
                end
                RRest=RRcmp+RRczest;
            end
        end
        QRSrefj=QRSrefj+round(RRest);
        j=j+1;     % increasing the index of qrs pointer vector
    else
        QRSrefj=QRSrefjlmi+imaxs+wleft;
    end
    
end
end %== function ================================================================
%
