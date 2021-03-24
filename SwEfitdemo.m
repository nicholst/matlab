% SwEdemo code
% Attempting to find minimum run time for mass-univariate application of sandwich estimator
%
% Based on https://github.com/nicholst/tenR/blob/master/SwEdemo.R
%
% T. Nichols 24 March 2021
% See https://github.com/nicholst/matlab/blob/master/LICENSE

% 'block' is cluster variable, might be family, site, subject (for repeated mes)
Nblock    = 15;
Nperblock = 4;       % Number of observations pers block
Nelm      = 32492;   % Number of vertices/voxels
rho       = 0.95;    % Intrablock correlation... maxed out to verify SwE is working


P       = 10;      % Number of predictors... all fake/simulated
N       = Nblock*Nperblock;

% Design: intercept, one between block variable, rest within block
X = [ones(N,1),...
     repelem(rand(Nblock,1),Nperblock),...
     rand(N,P-2)];

% Simulate repeated measures data, N x Nelm, with intrablock correlation rho
Y = sqrt(1-rho)*randn(N,Nelm) + ...
    sqrt(rho)  *reshape(repelem(randn(Nblock*Nelm,1),Nperblock),N,Nelm);

tic;

% OLS fit
pX      = pinv(X);
bh      = pX*Y;
res     = Y-X*bh;
sig2ols= sum(res.^2)/(N-P);
SEols  = diag(pX*pX') .* sig2ols;
Tols   = bh./SEols;

disp('OLS: ');toc


% Computation of SwE standard errors
tic;

Iblock = repelem([1:Nblock]',Nperblock);
[cbetahat,cbetaSE]=SwEfit(X,Iblock,Y);

disp('SwE vectorised: ');toc

Tsd=[std(Tols(2,:)),std(Tswe(2,:))];

disp('Standard deviation of T scores for between block covariate... should be 1.0\n')
fprintf('SD(T_ols) = %f  SD(T_sd) = %f\n',Tsd);





