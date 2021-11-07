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
Nelm      = 2000;    % Number of vertices/voxels
rho       = 0.95;    % Intrablock correlation... maxed out to verify SwE is working

P         = 10;      % Number of predictors... all fake/simulated
N         = Nblock*Nperblock;
Iblock    = repelem([1:Nblock]',Nperblock);
prop0     = 0.95;    % proportion of predictor that piles up at min value
alph      = 0.05;    % Nominal alpha for FPR report
nWB       = 100;    % Number of Wild Bootstrap iterations - Start at 100; try 1000

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

if nWB>0
    tic
    % Wild Bootstrap
    res0 = Y-X*cbetahat0;
    res1 = Y-X*cbetahat1;
    Pwb0 = zeros(P,Nelm); 
    Pwb1 = zeros(P,Nelm);
    for i=1:nWB
        if rem(i,10)==0; fprintf('%d ',i); end
        WBf=kron(2*binornd(1,0.5,Nblock,1)-1,ones(Nperblock,1));

        Ywb0 = WBf.*res0;
        [cbwb,cbSEwb] = SwEfit0(X,Iblock,Ywb0);
        Pwb0 = Pwb0 + (cbwb./cbSEwb >= Tswe0);

        Ywb0 = WBf.*res1;
        [cbwb,cbSEwb] = SwEfit(X,Iblock,Ywb0,[],1);
        Pwb1 = Pwb1 + (cbwb./cbSEwb >= Tswe1);
    end
    Pwb0 = (Pwb0+1)/(nWB+1);
    Pwb1 = (Pwb1+1)/(nWB+1);
    fprintf('\nSwE (both) Wild Bootstrap: ');toc
end    

cname={'Mean','PureBtwCov','PureWtnCov','BtwWtnCov'};

fprintf('\nEfficiency of non-iid working cov relative to OLS\n')
for i = 2:4
    fprintf('%12s: SD(beta_swe1)/SD(beta_ols) = %f\n',...
        cname{i},std(cbetahat1(i,:))/std(bh(i,:)));
end

fprintf('\nStandard deviation of T (should be 1.0)\n')
for i = 2:4
    fprintf('%12s: SD(T_ols) = %f  SD(T_swe0) = %f  SD(T_swe1) = %f\n',...
        cname{i},std(Tols(i,:)), std(Tswe0(i,:)),  std(Tswe1(i,:)));
end

Ta=tinv(1-alph,N-P);
fprintf('\nFPR (nominal %g, CI [%.4f,%.4f])\n',alph,alph+[-1,1]*sqrt(alph*(1-alph)/Nelm));
for i = 2:4
    fprintf('%12s: T_ols: %.3f  T_swe0: %.3f  T_swe0wb: %.3f  T_swe1: %.3f  T_swe1wb: %.3f\n',...
        cname{i},mean(Tols(i,:)>=Ta),mean(Tswe0(i,:)>=Ta),mean(Pwb0(i,:)<=alph),mean(Tswe1(i,:)>=Ta),mean(Pwb1(i,:)<=alph));
end

figure;imagesc(X)

Str={'Pure between','Pure within','Mixed'};
for i = 2:4
    figure
    P0 = tcdf(Tswe0(i,:),N-P,'upper');
    P1 = tcdf(Tswe1(i,:),N-P,'upper');
    loglog((1:Nelm)'/Nelm,[sort(P0'),sort(P1'),sort(Pwb0(i,:)'),sort(Pwb1(i,:)')],'linewidth',2);
    grid on;legend('swe0','swe1','sweWB0','sweWB1','AutoUpdate','off')
    set(refline(1),'linestyle',':')
    title(sprintf("Parameter %d: %s",i,Str{i-1}))
end


                                       
