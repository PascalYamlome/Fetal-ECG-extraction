function FiltECG= BandPassZeroPhasefilt(ECG, CF1,CF2,Filtorder,Sr)

    if size(ECG,1)>size(ECG,2)
        ECG = ECG';
    end
    
    bpfilter = designfilt('bandpassfir','FilterOrder',Filtorder, ...
         'CutoffFrequency1',CF1,'CutoffFrequency2',CF2, ...
         'SampleRate',Sr);
    
    for i = 1:size(ECG,1)-1
        FiltECG(i,:)= filtfilt(bpfilter,ECG(i,:));
    end
    FiltECG(6,:)=ECG(6,:);

end