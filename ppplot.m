function ppplot(Pval,varargin)
% FORMAT ppplot(Pval)
% FORMAT ppplot(Pval,<plot options>)
%
% Create -log10 PP-plot with pointwise confidence bands
%______________________________________________________________________________
% Author: T. Nichols
% Version: http://github.com/nicholst/matlab/tree/$Format:%h$
%          $Format:%ci$

% PP-Plot
N=length(Pval);
I=1:N;
plot(-log10(I/(N+1)),-log10(sort(Pval)),varargin{:})
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

h=get(gca,'children');
hold on;                                                                    
fill(-log10(xx),-log10(yy),[1 1 1]*0.85,'linestyle','none','facealpha',.5)
hold off                                                                    
uistack(h,'top');
title('-log10 P-P plot')

