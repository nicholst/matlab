%function SwEdemo(X,fid,Y,con)
% SwEdemo code
% Attempting to find minimum run time for mass-univariate application of sandwich estimator
%
% Based on https://github.com/nicholst/tenR/blob/master/SwEdemo.R
%
% T. Nichols 24 March 2021


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

disp("OLS: ");toc


% Computation of SwE standard errors
tic;

SEswe = zeros(P,Nelm);

Iblock = repelem([1:Nblock]',Nperblock);
Bread  = inv(X'*X)';
BreadX = Bread * X';
S      = zeros(P,P,Nelm);
S0     = zeros(1,P,Nelm);
for s = 1:Nblock
    I    = (s==Iblock);
    Ns   = sum(I);
    % BreadX times half of Meat
    S0   = BreadX(:,I)*res(I,:);
    S0   = reshape(S0,[1,P,Nelm]);
    % Full `Bread*Meat*Bread' contribution for block s
    S    = S + mtimesx(S0,'t',S0);
end

S = reshape(S,P*P,[]);
SEswe = sqrt(S(1:(P+1):end,:));
Tswe = bh./SEswe;

disp("SwE vectorised: ");toc

Tsd=[std(Tols(2,:)),std(Tswe(2,:))];

disp("Standard deviation of T scores for between block covariate... should be 1.0\n")
fprintf("SD(T_ols) = %f  SD(T_sd) = %f\n",Tsd);





