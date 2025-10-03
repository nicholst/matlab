function ppplot(Pval,varargin)
% FORMAT ppplot(Pval,[M],<plot options>)
%
% Create -log10 PP-plot with pointwise confidence bands
%
% M - logical vector of length(Pval), marking p-values to be highlighted; 
%     set to [] to pass plot options
%
%______________________________________________________________________________
% Author: T. Nichols
% Version: http://github.com/nicholst/matlab/tree/$Format:%h$
%          $Format:%ci$

if length(varargin)>0
    M=varargin{1};
    varargin(1)=[];
end

% Check for NaNs
if any(isnan(Pval))
    I=isnan(Pval);
    fprintf('WARNING: Removing NaN p-values (%d, %0.1f%%)\n',sum(I),mean(I)*100)
    Pval=Pval(~isnan(Pval));
    if length(Pval)==0
        error('All p-values NaN!')
    end
end

% PP-Plot
N=length(Pval);
I=1:N;
XX=-log10(I/(N+1));
YY=-log10(sort(Pval));
plot(XX,YY,varargin{:})
set(refline(1),'linestyle',':')
xlabel('Null Expected -log10 P')
ylabel('Observed -log10 P')

% Confidence Bands
Ncb=100;
alph=0.05;
bI=10.^[linspace(log10(1),log10(N),Ncb)'];
rI=flipud(bI);
xx=[bI;rI]/(N+1);
% https://en.wikipedia.org/wiki/Order_statistic#Order_statistics_sampled_from_a_uniform_distribution
yy=[betainv(alph/2,bI,N-bI+1);betainv(1-alph/2,rI,N-rI+1)];

if ~isempty(M)
    if all(size(Pval)~=size(M))
        error("Pvalue and mask vector don't match")
    end
    M=logical(M);
    [~,I]=sort(Pval);
    line(XX(M(I)),YY(M(I)),'marker','o','linestyle','none');
end

h=get(gca,'children');
hold on;                                                                    
fill(-log10(xx),-log10(yy),[1 1 1]*0.85,'linestyle','none','facealpha',.5)
hold off                                                                    
uistack(h,'top');
title('-log10 P-P plot')

