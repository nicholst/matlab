function [pID,pN,pA] = FDR(p,aFDR)
% FORMAT [pID,pN,pA] = FDR(p,aFDR)
% 
% p     - vector of p-values
% aFDR  - alpha_FDR, aka q, the desired False Discovery Rate level
%
% pID - p-value threshold based on independence or positive dependence
% pN  - Nonparametric p-value threshold
% pA  - Adaptive FDR result of Benjamini, Krieger, Yekutieli, 2006. Adaptive 
%       linear step-up procedures that control the false discovery rate. 
%       Biometrika 93, 491–507. doi:10.1093/biomet/93. 
%______________________________________________________________________________
% Author: T. Nichols
% Version: http://github.com/nicholst/matlab/tree/$Format:%h$
%          $Format:%ci$


p = p(isfinite(p));  % Toss NaN's
p = sort(p(:));
V = length(p);
I = (1:V)';

cVID = 1;
cVN = sum(1./(1:V));
denomVA = V + 1 - I*(1-aFDR);


pID = p(find(p<=I/V*aFDR/cVID,1,'last'));
if isempty(pID)
  pID=0;
end

pN = p(find(p<=I/V*aFDR/cVN,1,'last'));
if isempty(pN)
  pN=0; 
end

% This is non-looping implimentation of the BKY method... It certainly
% doesn't reject too much, but need to double-check I'm not missing any
% rejections.
tmp = find(p>I*aFDR/denomVA,1,'first');
tmp = max(tmp-1,1);
pA = p(tmp);
if isempty(pA)
  pA=0;
end
