function x = meansc(v,perci,percf)
% --------------------------------------------------------------------------------------------
% meansc.m: Compute the mean value of a vector excluding the distribution tails.
%   x= meansc(v,perci,percf)
%           "perci" = % of min values;  "percf" = % of max values to be excluded
%
%   Version 1.00, Date: 08/04/2000
% --------------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% --------------------------------------------------------------------------------------------

if(nargin<2), perci=5; end
if(nargin<3), percf=perci; end
if(perci==50  && percf==50), x=median(v); return; end
if(perci+percf>=100), x=[]; return; end
vo=sort(v);
ii = 1+floor(length(v)*perci/100);
fi = length(v)-floor(length(v)*percf/100);
x=mean(vo(ii:fi));

end %== function ================================================================
%
