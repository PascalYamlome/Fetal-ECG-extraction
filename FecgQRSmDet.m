function qrsM=FecgQRSmDet(Se,fs,cName,graph,dbFlag,saveFig,qrsAf)
% ---------------------------------------------------------------------------------------------
%   Fecg: "Mother" QRS detection
%  - Best mother ECG selection based on a priori information on mother ECG pseudo-periodicity.
%    A raw derivative filter signal obtained as the difference between 
%    the average values on two intervals of 7ms far off 9ms)
%  - The absolute value of the raw derivative signal was filtered by a forward backward 
%    Butterworth bandpass filter (6.3-16.Hz).
%  - mother QRS detection based on this absolute derivative
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
% --------------------------------------------------------------------------------------------

if(nargin<3), cName=''; end
if(nargin<4), graph=1; end
if(nargin<5), dbFlag=0; end
if(nargin<6), saveFig=0; end
if(nargin<7), qrsAf=[]; end

fprintf('\n --------------------------------------------------------- \n');
[progpath, progname] = fileparts(which(mfilename));
fprintf('Program: %s,   record name: %s\n', progname, cName);

if(isempty(qrsAf));
    learning=0;   % test set, fetal annotations are not available
else
    learning=1;   % learning set, fetal annotations are available
end

%-------------------------------------------------------------
% recording duration
[ndt, ns]=size(Se);
vtime= [1:ndt]/fs;

if(dbFlag && learning), PlotSgnMrkNc(X, qrsAf*fs, fs, cName); end

for is=1:ns,
    X(:,is)= (Se(:,is)-mean(Se(:,is)))/std(Se(:,is));
end
if(graph)
    figure, set(gcf,'Color','white');
    for is=1:ns,
        subplot(ns,1,is), plot(vtime,X(:,is));
        wgmi1= min(X(:,is)) -2;
        wgma1= max(X(:,is)) +2;
        ylim([wgmi1, wgma1]);
        set(gca,'YTick',[-5 0 5])
        % if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': ECG signals']); end
    end
end
% ----
ecg=X;
ecgf=ecg;

% raw derivative filter coefficients (7ms before and after, 9ms in the
% middle)
nu=ceil(0.0070 * fs); nz=floor(0.0090*fs /2)*2 +1;  % nz=nearest odd value
B=[ones(nu,1);zeros(nz,1);-ones(nu,1)];
delay=floor(length(B)/2);
% compute the derivative signal
ecgfx=[repmat(ecgf(1,:),delay,1);ecgf;repmat(ecgf(end,:),delay,1)];
decgr=filter(B,1,ecgfx);   decgr= decgr(2*delay+1:end,:);

adecg=abs(decgr); % absolute ECG derivative

% -----------    chose the best mother ecg channel
w8=fix(8*fs);              % window of 8 s
w2=fix(2*fs);              % window of 2 s
w02=fix(0.2*fs);           % window of 0.2 s
mD8=zeros(ns,1); mD2=zeros(ns,1); mD02=zeros(ns,1);
for is=1:ns, 
    % compute the average of maximum derivatives on windows of 8s (5% of maxima are discarted)
    mD8(is)=meanMaxSc(adecg(:,is), w8, 0,5);
    % compute the average of maximum derivatives on windows of 2s (5% of maxima are discarted)
    mD2(is)=meanMaxSc(adecg(:,is), w2, 0,5);  
    % compute the average of maximum derivatives on windows of 0.2s (1% of maxima are discarted)
    mD02(is)=meanMaxSc(adecg(:,is), w02, 0,1);  
end
qualFact=mD2./(mD02+mD8);
[dummy,ics]=sort(-qualFact);
% ------------
Ses=Se(:,ics);
ecgs=ecg(:,ics);
decgs=decgr(:,ics);
adecgs=adecg(:,ics);
% ----
if(graph || saveFig)
    figure, set(gcf,'Color','white');
    for is=1:ns,
        subplot(ns,1,is), plot(vtime, Ses(:,is));
        wgmi1= min(Ses(:,is)) -2;
        wgma1= max(Ses(:,is)) +2;
        ylim([wgmi1, wgma1]);
        set(gca,'YTick',[-5 0 5])
        % if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': sorted mother ECG']); end
        shg
    end
    if(saveFig), figFmt='png';
        figPath=fullfile('../Figure/',progname);
        if(~exist(figPath,'dir')), mkdir(figPath); end
        figName=fullfile(figPath,[cName,'_mother']);
        print(gcf, ['-d',figFmt],figName);
    end
end

if(dbFlag && learning), PlotSgnMrkNc(Ses', qrsAf*fs, fs, [cName,' - eS']); end
% ----
ecg1=ecgs(:,1);
decg1=decgs(:,1);
adecg1=adecgs(:,1);

% The absolute value of the raw derivative signal was filtered by 
% a forward-backward  Butterworth bandpass filter (6.3-16.Hz).
fmind=5; fmaxd=20;
Wn = [fmind, fmaxd]/(fs/2);  % normalized cut-off frequency
[b,a]= butter(1,Wn);
adecg1=filtfilt(b,a,adecg1);


npx=fix(1*fs);
decg1x=[zeros(npx,1);decg1;zeros(npx,1)];
adecg1x=[zeros(npx,1);adecg1;zeros(npx,1)];

% --- QRS detection -----
% qrsM(1,:) QRSref    (reference point, max signed derivative )
% qrsM(2,:) QRSonset  (derivative overcomes half threshold)
% qrsM(3,:) supThDer  (derivative overcomes threshold)
% qrsM(4,:) maxAbsDer (max absolute derivative)
% qrsM(5,:) infThDer  (derivative becomes lower than threshold)
% qrsM(5,:) QRSoffset (derivative decreses below half threshold)
pth=0.45; qrsM=QRSdetectorM(adecg1x,decg1x,fs,pth,0.85,1);

qrsM=qrsM-npx;
qrsRef=qrsM(1,:);
nQRS=length(qrsRef);
fprintf('Number of detected mother QRSs= %d\n', nQRS);

RRc= diff(qrsRef);  RRs= RRc/fs;

if(learning), RRafs= diff(qrsAf); end    % annotated fetal qrs
[nb,vp]=hist(RRs,[0.4:0.05:1.9]);
if(graph)
    figure, bar(vp,nb,1); title([cName,': detected RR hist']);
end
% [nbm,ibm]=max(nb);  RRst=vp(ibm);   % most frequent value
% iqrsiT=find(RRst-0.05<RRs & RRs<RRst+0.05);

RRmean=meansc(RRs,4,4);
RRstd= std(RRs);
fprintf('RR mean= %5.3f,  stdev=%5.3f\n', RRmean, RRstd);
if(graph && dbFlag)
    figure
    subplot(2,1,1),plot(adecg1,'b');
    %    DrawVertMarker(qrsM(2,:)','c',':','none');
    %    DrawVertMarker(qrsM(end,:)','c',':','none');
    DrawVertMarker(qrsM(1,:)','r',':','none');           % reference point
    subplot(2,1,2),plot(ecg1,'b');
    %    DrawVertMarker(qrsM(2,:)','c',':','none');
    %    DrawVertMarker(qrsM(end,:)','c',':','none');
    DrawVertMarker(qrsM(1,:)','r',':','none');           % reference point
end

if(graph || saveFig)
    figure, set(gcf,'Color','white'); hold on;
    if(learning), plot(qrsAf, [RRafs(1);RRafs],'r.-'); end
    plot(qrsM(1,:)/fs, [RRs(1),RRs],'b.-');
    ylim([0.3, 2]);
    if(learning), title([cName,': mother (b) & annotated fetal(r) RR series']);
    else title([cName,': mother (b) RR series']); end
    figResize(0, 1, 1, .35);
    if(saveFig), figFmt='png';
        figPath=fullfile('../Figure/',progname);
        if(~exist(figPath,'dir')), mkdir(figPath); end
        figName=fullfile(figPath,[cName,'_RRmother']);
        print(gcf, ['-d',figFmt],figName);
    end
end
if(graph)
    figure, set(gcf,'Color','white');
    for is=1:ns,
        subplot(ns,1,is), plot(vtime, Ses(:,is));
        wgmi1= min(Ses(:,is)) -2;
        wgma1= max(Ses(:,is)) +2;
        ylim([wgmi1, wgma1]);
        %   DrawVertLines(vtime', qrsM(1,:)/fs, 'r', ylim);
        DrawVertMarker(qrsM(1,:)'/fs,'r',':','none');
        %     DrawVertMarker(qrsM(2,:)'/fs,'b',':','none');
        %     DrawVertMarker(qrsM(end,:)'/fs,'r',':','none');
        set(gca,'YTick',[-5 0 5])
        % if(is~=ns), set(gca,'XTickLabel',''); end
        if(is==1), title([cName,': ECGs with mother QRS marker ']); end
        shg
    end
end
if(graph || saveFig)
    PlotSgnMrkNc(ecgs(:,1), qrsM(1,:), fs, cName);
    shg
    if(saveFig), figFmt='png';
        figPath=fullfile('../Figure/',progname);
        if(~exist(figPath,'dir')), mkdir(figPath); end
        figName=fullfile(figPath,[cName,'_mother']);
        print(gcf, ['-d',figFmt],figName);
    end
end

return
end %== function ================================================================
%

