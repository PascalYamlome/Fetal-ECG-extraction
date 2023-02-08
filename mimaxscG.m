% --------------------------------------------------------------------------------------------
% mimaxscG.m: Compute the minimum and the maximum value of a vector excluding
%             the distribution tails.
%   mimaV = mimaxscG(v,perci,percf,marg)
%   "perci" = % of min values to be excluded  
%   "percf" = % of max values to be excluded
%   "marg"  = margin as fraction of (max-min), subctracted to mimaV(1) and added to mimaV(2) 
%   "mimaV" = output two element vector containing the min and the max
%
% --------------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% --------------------------------------------------------------------------------------------
function mimaV = mimaxscG(v,perci,percf,pmarg)

vo=sort(v(:));
if(nargin<2), perci=5; end
if(nargin<3), percf=perci; end
if(nargin<4), pmarg=0.1; end

if(perci<0), ii=1-perci;
else ii=1+floor(length(vo)*perci/100); end
if(percf<0), fi=length(vo)+percf;
else fi=length(vo)-floor(length(vo)*percf/100); end

mi=min(vo(ii:fi));
ma=max(vo(ii:fi));
marg=pmarg *(ma-mi);

mimaV=[mi-marg, ma+marg];

end %== function ================================================================
%
