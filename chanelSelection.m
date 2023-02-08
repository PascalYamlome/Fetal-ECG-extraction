function sel = chanelSelection(Se, fs)
    [ndt, ns]=size(Se);
    vtime= [1:ndt]/fs;
    
    %if(dbFlag && learning), PlotSgnMrkNc(X, qrsAf*fs, fs, cName); end
    %mean over sd for each channel scaled and mean normalised
    for is=1:ns,
        X(:,is)= (Se(:,is)-mean(Se(:,is)))/std(Se(:,is));
    end
    
%     figure, set(gcf,'Color','white');
%     for is=1:ns,
%         subplot(ns,1,is), plot(vtime,X(:,is));
%         wgmi1= min(X(:,is)) -2;
%         wgma1= max(X(:,is)) +2;
%         ylim([wgmi1, wgma1]);
%         set(gca,'YTick',[-5 0 5])
%         % if(is~=ns), set(gca,'XTickLabel',''); end
%         if(is==1), title([cName,': ECG signals']); end
%     end
%     
%  

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
pth=0.45; qrsM=QRSdetectorM(adecg1x,decg1x,fs,pth,0.85,1);

qrsM=qrsM-npx;
qrsRef=qrsM(1,:);
nQRS=length(qrsRef);
fprintf('Number of detected mother QRSs= %d\n', nQRS);

RRc= diff(qrsRef);  RRs= RRc/fs;

% if(learning), RRafs= diff(qrsAf); end    % annotated fetal qrs
% [nb,vp]=hist(RRs,[0.4:0.05:1.9]);
% if(graph)
%     figure, bar(vp,nb,1); title([cName,': detected RR hist']);
% end



    
    
    sel = qrsM;

end