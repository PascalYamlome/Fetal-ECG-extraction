function [mi, ma] = mimaxsc(v,perci,percf)
% --------------------------------------------------------------------------------------------
% mimaxsc.m: Compute the minimum and the maximum value of a vector excluding 
%            the distribution tails.
%   [mi, ma] = mimaxsc(v,perci,percf)
%    "v"     = input data vector
%    "perci" = % of min values to be excluded
%    "percf" = % of max values to be excluded
%
% --------------------------------------------------------------------------------------------
%   Maurizio Varanini, Clinical Physiology Institute, CNR, Pisa, Italy
%   For any comment or bug report, please send e-mail to: maurizio.varanini@ifc.cnr.it
% --------------------------------------------------------------------------------------------

vo=sort(v(:));
if(nargin<2), perci=5; end
if(nargin<3), percf=perci; end

if(perci<0), ii=1-perci;
else ii=1+floor(length(v)*perci/100); end
if(percf<0), fi=length(v)+percf;
else fi=length(v)-floor(length(v)*percf/100); end

mi=min(vo(ii:fi));
ma=max(vo(ii:fi));
end %== function ================================================================
%
