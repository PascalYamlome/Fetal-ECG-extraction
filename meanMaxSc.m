function me = meanMaxSc(v,nel,percmi,percma)
% --------------------------------------------------------------------------------------------
% meanMaxSc.m: Compute the average value of the maxima computed on data windows
%              of the input vector.
%  mv = meanMaxSc(v,nel,percmi,percma)
%              Distribution tails can be excluded using parameters "perci" and "percf"
%              "v" = input data vector
%              "nel" = data window length
%              "percmi" = % of min values to be excluded; (if positive)
%                         number of min values (if negative)
%              "percma" = % of max values to be excluded (if positive)
%                         number of max values (if negative)
%
%   Version 1.00, Date: 01/05/2000
% ---------------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% ---------------------------------------------------------------------------------------------

if(nargin<3), percmi=5; end
if(nargin<4), percma=percmi; end

if(nel>=length(v)), me=max(v); return; end

j=1;
maxi=zeros(floor(length(v)/nel),1);
for i=1:nel:length(v)-nel+1
    maxi(j)=max(v(i:i+nel-1));
    j=j+1;
end

if(percmi<0), ii=1-percmi;
else ii=1+floor(length(maxi)*percmi/100); end
if(percma<0), fi=length(maxi)+percma;
else fi=length(maxi)-floor(length(maxi)*percma/100); end

omaxi=sort(maxi);
me=mean(omaxi(ii:fi));

end %== function ================================================================
%
