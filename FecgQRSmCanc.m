function Xr=FecgQRSmCanc(X,qrsM,fs,cName);%,graph,dbFlag,saveFig,qrsAfs)
% ---------------------------------------------------------------------------------------------
% Fecg: "Mother" ECG canceling using  PQRST approximation obtained by 
% Singular Value Decomposition (SVD).
% A trapezoidal window is used to select and weight the signal around each detected mother QRS.
%
% Author: Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
% For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% ---------------------------------------------------------------------------------------------
% This program is free software; you can redistribute it and/or modify it under the terms
% of the GNU General Public License as published by the Free Software Foundation; either
% version 2 of the License, or (at your option) any later version.
%
% This program is distributed "as is" and "as available" in the hope that it will be useful,
% but WITHOUT ANY WARRANTY of any kind; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% ---------------------------------------------------------------------------------------------

if(nargin<4), cName=''; end
if(nargin<5), graph=1; end
if(nargin<6), dbFlag=1; end
if(nargin<7), saveFig=0; end
if(nargin<8), qrsAfs=[]; end

grafEigD=0;

fprintf('\n --------------------------------------------------------- \n');
[progpath, progname] = fileparts(which(mfilename));
fprintf('Program:  %s,  record name: %s\n', progname, cName);

if(~isempty(qrsAfs));
    learning=1;   % learning set, annotations are available
else
    learning=0;   % test set
end

%-------------------------------------------------------------
% recording duration
[ndt, ns]=size(X);

if(dbFlag && learning), PlotSgnMrkNc(X, qrsAfs*fs, fs, [cName, ' - fetal QRS ann.']); end

vtime= [1:ndt]/fs;
% -----
% qrsM(1,:) QRSref    (reference point, max signed derivative )
% qrsM(2,:) QRSonset  (derivative overcomes half threshold)
% qrsM(3,:) supThDer  (derivative overcomes threshold)
% qrsM(4,:) maxAbsDer (max absolute derivative)
% qrsM(5,:) infThDer  (derivative becomes lower than threshold)
% qrsM(5,:) QRSoffset (derivative decreses below half threshold)

qrsR=qrsM(1,:);        % the reference point is the point of max signed derivative
% if(dbFlag), PlotSgnMrkNc(X, qrsR, fs, [cName, ' - mother QRS det.']); end
nQRS=length(qrsR);
fprintf('Number of QRSs= %d\n', nQRS);

RRc= diff(qrsR);
RRs= RRc/fs;
RRmean=meansc(RRs,4,4);
%RRmean= mean(RRs);
RRstd= std(RRs);
fprintf('RR mean= %d,  stdev=%d\n', RRmean, RRstd);

% the reference point is the point of max signed derivative
npp=fix(0.2*fs);               % number of samples before the qrs reference
npd=fix(min(0.5, 0.8*(RRmean-0.1))*fs);  % number of samples after the qrs reference
npt=1+npp+npd;

% extend the signals in order to manage the first and the last QRS
Xx=X;
npqp=fix(0.12*fs);
if(qrsR(1)-npqp < 1), qi=2; 
    npxp=0;
    Xx(1:npqp,:)=repmat(Xx(npqp+1,:),npqp,1);
else qi=1;
    npxp=max(0,npp+1-qrsR(1));  % number of samples added to left of each signal
    Xx=[repmat(X(1,:),npxp,1);X];
end
npqd=fix(0.14*fs);
qf=nQRS;
if(qrsR(end)+npqd > size(X,1)), qf=qf-1;
    Xx(end-npqd+1:end,:)=repmat(Xx(end-npqd-1,:),npqd,1);
end
if(qrsR(qf)+0.85*fs*median(RRs(end-4:end)) < size(X,1))
    npep=fix(max(0.1*fs, 0.15*fs*mean(RRs(end-4:end))));
    Xx(end-npep+1:end,:)=repmat(Xx(end-npep-1,:),npep,1);
end
npxd=max(0,qrsR(qf)+npd-ndt);  % number of samples added to right of each signal
Xx=[Xx;repmat(X(end,:),npxd,1)];

ndtx=size(Xx,1);
nqe= qf-qi+1;
vtimex= [1-npxp:ndt+npxd]/fs;


% built vectors of indexes (start and end) of QRS window
iqw=qrsR(qi:qf);
iw=npxp+iqw-npp;      %  start of QRS window
fw=npxp+iqw+npd;      %  end of QRS window
A=zeros(npt,nqe);
wwg=weightFun2(npp,npd,fs);
Xc=zeros(ndtx,ns);
for is=1:ns
    for iq=1:nqe
        iwq=iw(iq); fwq=fw(iq);
        A(:,iq)=Xx(iwq:fwq,is).*wwg;
    end
    % Singular value decomposition
    % S = diagonal matrix containing the singular values
    % U and V = unitary matrices so that X = U*S*V'.
%    tic;
    [U,S,V] = svd(A,0);
%    fprintf('SVD (%d,%d) computation time = %3.2f\n', size(A), toc);
    
%    fprintf('%f,  ', diag(S));   fprintf('\n');
    if(grafEigD), figure; plot(diag(S),'r-*');
        title([num2str(is),'- Singular values']); drawnow;
    end
    
    mt=nqe;
        
    nds=3;     % select the number of singular values  <===
    
    sv=diag(S);
    pownds(is)=sum(sv(1:nds).^2)/sum(sv(1:nqe).^2);
    if(dbFlag), fprintf('is=%d, pownds=%f\n', is, pownds(is)); end
    
%    tic
    Sr=S; for k=1:mt-nds-1, Sr(mt-k,mt-k)=0; end
    Ar=U*Sr*V';
%    fprintf('Ar (%d,%d) rebuilding time = %3.2f\n', size(Ar), toc);
    
    for iq=1:nqe
        iwq=iw(iq); fwq=fw(iq);
%        Xc(iwq:fwq,is)=Ar(:,iq);
        Xc(iwq:fwq,is)=Ar(:,iq)./wwg;   % .*windLR;
    end
    
    % smoothing connections between successive windows
    for iq=1:nqe-1
        fwq=fw(iq); iwqs=iw(iq+1);
        if(iwqs>fwq)
            dv=Xc(iwqs,is)-Xc(fwq,is);
            pv=dv/(iwqs-fwq);
            Xc(fwq+1:iwqs-1,is)=Xc(fwq,is)+pv*(1:iwqs-fwq-1)';
        end
    end
end
% if(dbFlag)
%     for is=1:ns,
%         figure
%         set(gcf,'Color','white');
%         subplot(2,1,1), plot(vtimex,[Xx(:,is),Xc(:,is)]);
%         wgmi1= min(Xx(:,is)) -2;
%         wgma1= max(Xx(:,is)) +2;
%         ylim([wgmi1, wgma1]);
%         set(gca,'YTick',[-5 0 5])
%         DrawVertMarker(iw'/fs,'r',':','none');
%         DrawVertMarker(qrsM(1,:)'/fs,'g',':','none');
%         DrawVertMarker(fw'/fs,'m',':','none');
%         title([cName,' - ',num2str(is),': ecg & ecgM']);
%         subplot(2,1,2), plot(vtimex,[Xx(:,is)-Xc(:,is)]);
%         ylim([wgmi1, wgma1]);
%         set(gca,'YTick',[-5 0 5])
%         DrawVertMarker(iw'/fs,'r',':','none');
%         DrawVertMarker(qrsM(1,:)'/fs,'g',':','none');
%         DrawVertMarker(fw'/fs,'m',':','none');
%         title([cName,' - ',num2str(is),': ecg - ecgM']);
%         shg
%     end
% end
Xc= Xc(npxp+1:end-npxd,:);
Xx= Xx(npxp+1:end-npxd,:);
Xe=Xx-Xc;
iw=iw-npxp;  fw=fw-npxp;
% if(graph)
%     for is=1:ns,
%         figure,  set(gcf,'Color','white');
%         subplot(2,1,1), plot(vtime,[X(:,is),Xc(:,is)]);
%         wgmi1= min(X(:,is)) -2;
%         wgma1= max(X(:,is)) +2;
%         ylim([wgmi1, wgma1]);
%         set(gca,'YTick',[-5 0 5])
%         DrawVertMarker(iw'/fs,'r',':','none');
%         DrawVertMarker(qrsM(1,:)'/fs,'g',':','none');
%         DrawVertMarker(fw'/fs,'m',':','none');
%         title([cName,' - ',num2str(is),': ecg & ecgM']);
%         subplot(2,1,2), plot(vtime,Xe(:,is));
%         ylim([wgmi1, wgma1]);
%         set(gca,'YTick',[-5 0 5])
%         DrawVertMarker(iw'/fs,'r',':','none');
%         DrawVertMarker(qrsM(1,:)'/fs,'g',':','none');
%         DrawVertMarker(fw'/fs,'m',':','none');
%         title([cName,' - ',num2str(is),': ecg - ecgM']);
%         shg
%     end
% end
% %--------------------------------------------------------------------------------------------
% if(graph || saveFig)
%     figure, set(gcf,'Color','white');
%     for is=1:ns,
%         subplot(ns,1,is), plot(vtime,Xe(:,is));
%         wgmi1= min(Xe(:,is)) -2;
%         wgma1= max(Xe(:,is)) +2;
%         ylim([wgmi1, wgma1]);
%         DrawVertMarker(iw'/fs,'r',':','none');
%         DrawVertMarker(qrsM(1,:)'/fs,'g',':','none');
%         DrawVertMarker(fw'/fs,'m',':','none');
%         set(gca,'YTick',[-5 0 5])
%         % if(is~=ns), set(gca,'XTickLabel',''); end
%         if(is==1), title([cName,': mixed foetal ECG']); end
%     end
%     shg
%     if(saveFig), figFmt='png';
%         figPath=fullfile('../Figure/',progname);
%         if(~exist(figPath,'dir')), mkdir(figPath); end
%         figName=fullfile(figPath,[cName,'_mCanc']);
%         print(gcf, ['-d',figFmt],figName);
%         %    saveas(gcf, figName,'fig');
%     end
% end
%-------------------------------------------------------------

Xr=Xe;

return
end % = function ===============================================================
% ------------------------------------------------------------------
function wwg=weightFun2(npp,npd,fs)
nppc1=fix(0.06*fs); npdc1=fix(0.06*fs);
nppc2=fix(0.08*fs); npdc2=min(fix(0.2*fs),npd-npdc1);
ii1=1; ie1=npp-nppc1-nppc2;
ii2=ie1+1; ie2=ie1+nppc2;
ii3=ie2+1; ie3=ie2+nppc1+npdc1+1;
ii4=ie3+1; ie4=ie3+npdc2;
ii5=ie4+1; ie5=npp+npd+1;
wwg(ii1:ie1)=0.20;
wwg(ii2:ie2)=0.20+0.8*(1:nppc2)/nppc2;
wwg(ii3:ie3)=1;
wwg(ii4:ie4)=1-0.8*(1:npdc2)/npdc2;
wwg(ii5:ie5)=0.20;
wwg=wwg';
% figure, plot(wwg);
return
end % = function ===============================================================
%

