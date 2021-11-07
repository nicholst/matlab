% SwEdemo code
% Attempting to find minimum run time for mass-univariate application of sandwich estimator
%
% Based on https://github.com/nicholst/tenR/blob/master/SwEdemo.R
%
% T. Nichols 27 March 2021
% See https://github.com/nicholst/matlab/blob/master/LICENSE

% 'block' is cluster variable, might be family, site, subject (for repeated mes)
Nblock    = 1000;
Nperblock = 4;       % Number of observations pers block
Nelm      = 10000;    % Number of vertices/voxels
rho       = 0.95;    % Intrablock correlation... maxed out to verify SwE is working

P         = 10;      % Number of predictors... all fake/simulated
N         = Nblock*Nperblock;
Iblock    = repelem([1:Nblock]',Nperblock);
prop0     = 0.95;     % proportion of predictor that piles up at min value
alph      = 0.001;    % Nominal alpha for FPR report

% Design: intercept, one between block variable, rest within block
myrand=@(p,q,prop0)[rand(p,q).*binornd(1,1-prop0,p,q)];
%myrand=@(p,q,prop0)[(0.5+0.5*rand(p,q)).*binornd(1,1-prop0,p,q)];
X = [ones(N,1),...
     zeros(N,1),...
     myrand(N,P-2,prop0)];


% Split predictor 3: 
%     X(:,2) is pure between block
%     X(:,3) is pure within block
for i=1:Nblock
    I = Iblock==i;
    X(I,2) = mean(X(I,3));
    X(I,3) = X(I,3)-X(I,2);
    if binornd(1,prop0)
        X(I,2) = -0.5;
    end
end

% Simulate repeated measures data, N x Nelm, with intrablock correlation rho
Y = sqrt(1-rho)*randn(N,Nelm) + ...
    sqrt(rho)  *reshape(repelem(randn(Nblock*Nelm,1),Nperblock),N,Nelm);

% OLS fit
tic;
pX      = pinv(X);
bh      = pX*Y;
res     = Y-X*bh;
sig2ols = sum(res.^2)/(N-P);
SEols   = sqrt(diag(pX*pX') .* sig2ols);
Tols    = bh./SEols;
fprintf('OLS:                      ');toc


% Computation of SwE standard errors, iid working cov
tic;
[cbetahat0,cbetaSE0]=SwEfit0(X,Iblock,Y);
fprintf('SwE vectorised, ident W:  ');toc

% Computation of SwE standard errors, global working cov
tic;
[cbetahat1,cbetaSE1,Vg]=SwEfit(X,Iblock,Y,[],1);
fprintf('SwE vectorised, global W: ');toc

Tswe0 = cbetahat0./cbetaSE0;
Tswe1 = cbetahat1./cbetaSE1;

fprintf('\nEfficiency of non-iid working cov relative to OLS\n')
fprintf('PureBtwCov: SD(beta_swe1)/SD(beta_ols) = %f\n',...
        std(cbetahat1(2,:))/std(bh(2,:)));
fprintf('PureWtnCov: SD(beta_swe1)/SD(beta_ols) = %f\n',...
        std(cbetahat1(3,:))/std(bh(3,:)));
fprintf('BtwWtnCov:  SD(beta_swe1)/SD(beta_ols) = %f\n',...
        std(cbetahat1(4,:))/std(bh(4,:)));

fprintf('\nStandard deviation of T (should be 1.0)\n')
fprintf('PureBtwCov: SD(T_ols) = %f  SD(T_swe0) = %f  SD(T_swe1) = %f\n',...
        std(Tols(2,:)), std(Tswe0(2,:)),  std(Tswe1(2,:)));
fprintf('PureWtnCov: SD(T_ols) = %f  SD(T_swe0) = %f  SD(T_swe1) = %f\n',...
        std(Tols(3,:)), std(Tswe0(3,:)),  std(Tswe1(3,:)));
fprintf('BtwWtnCov:  SD(T_ols) = %f  SD(T_swe0) = %f  SD(T_swe1) = %f\n',...
        std(Tols(4,:)), std(Tswe0(4,:)),  std(Tswe1(4,:)));

Ta=tinv(1-alph,N-P);
fprintf('\nFPR (nominal %g, CI [%.4f,%.4f])\n',alph,alph+[-1,1]*sqrt(alph*(1-alph)/Nelm));
fprintf('PureBtwCov: FPR(T_ols) = %f  FPR(T_swe0) = %f  FPR(T_swe1) = %f\n',...
        mean(Tols(2,:)>=Ta), mean(Tswe0(2,:)>=Ta),  mean(Tswe1(2,:)>=Ta));
fprintf('PureWtnCov: FPR(T_ols) = %f  FPR(T_swe0) = %f  FPR(T_swe1) = %f\n',...
        mean(Tols(3,:)>=Ta), mean(Tswe0(3,:)>=Ta),  mean(Tswe1(3,:)>=Ta));
fprintf('BtwWtnCov:  FPR(T_ols) = %f  FPR(T_swe0) = %f  FPR(T_swe1) = %f\n',...
        mean(Tols(4,:)>=Ta), mean(Tswe0(4,:)>=Ta),  mean(Tswe1(4,:)>=Ta));

Str={'Pure between','Pure within','Mixed'};
for i = 2:4
    figure(i)
    P0 = tcdf(Tswe0(i,:),N-P,'upper');
    P1 = tcdf(Tswe1(i,:),N-P,'upper');
    loglog((1:Nelm)'/Nelm,[sort(P0'),sort(P1')]);
    grid on;legend('swe0','swe1','AutoUpdate','off')
    set(refline(1),'linestyle',':')
    title(sprintf("Parameter %d: %s",i,Str{i-1}))
end

figure(5);nhist({Tswe0(4,:),Tswe1(4,:)});title('Parameter 4');legend({'swe0','swe1'})

