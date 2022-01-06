% SwEdemo code
% Attempting to find minimum run time for mass-univariate application of sandwich estimator
%
% Based on https://github.com/nicholst/tenR/blob/master/SwEdemo.R
%
% T. Nichols 27 March 2021
% See https://github.com/nicholst/matlab/blob/master/LICENSE

% 'block' is cluster variable, might be family, site, subject (for repeated mes)
Nelm      = 1000;    % Number of vertices/voxels
rho       = 0.95;  % Intrablock correlation... maxed out to verify SwE is working
alph      = 0.05;  % Nominal alpha for FPR report
nWB       = 0;   % Number of Wild Bootstrap iterations - Start at 100; try 1000

% Design: Actual ABCD 2-visit design...
%
%    12525 measurements
%           8039 baseline measurements
%           4486 year 2 measurements
%
%    8740 subjects
%           4254 subjects with only baseline
%           3785 subjects with both visits
%            701 subjects with only year 2
%
%    7493 familes
%       Family sizes (by scan)
%              1    2    3    4    5    6    7 
%           3675 3051  350  400    5   11    1 
%       Family sizes (by subject)
%              1    2    3    4 
%           6278 1184   30    1  
%    
X     = load('FamSubVisX.dat');
Fam   = X(:,1);
Sub   = X(:,2);
Vis   = X(:,3);
Age   = X(:,4);
X     = X(:,5:end);

N     = size(X,1);
P     = size(X,2);

Nms   = readcell('FamSubVisX_names.txt', 'Delimiter',"");
Nms(1:4)=[];

Iblock = Fam;
Nblock = length(unique(Fam));
Nperblock = groupcounts(Iblock);

% Need rows in family order to simplify simulation
[~,Oblock]=sort(Iblock);
X      = X(Oblock,:);
Fam    = Fam(Oblock,:);
Sub    = Sub(Oblock,:);
Vis    = Vis(Oblock,:);
Age    = Age(Oblock,:);

% Simulate repeated measures data, N x Nelm, with intrablock correlation rho
Y = sqrt(1-rho)*randn(N,Nelm) + ...
    sqrt(rho)  *repelem(randn(Nblock,Nelm),Nperblock,1);

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
        WBf=repelem(2*binornd(1,0.5,Nblock,1)-1,Nperblock,1);

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
else
    Pwb0 = repelem(NaN,P,Nelm);
    Pwb1 = repelem(NaN,P,Nelm);
end    

%cname={'Mean','PureBtwCov','BtwWtnCov'};  %,'PureWtnCov'
cname={'Mean','BaseAge','DeltaAge'};  %,'PureWtnCov'
cI=[size(X,2),           1,           2]; 

fprintf('\nEfficiency of non-iid working cov relative to OLS\n')
fprintf('               SD(beta_swe1)/SD(beta_ols)\n');
for i = 1:length(cI)
    c=cI(i);
    fprintf('%12s:             %.3f\n',...
        cname{i},std(cbetahat1(c,:))/std(bh(c,:)));
end

fprintf('\nStandard deviation of T (should be 1.0)\n')
fprintf('                SD(T_ols) SD(T_swe0) SD(T_swe1)\n');
for i = 1:length(cI)
    c=cI(i);
    fprintf('%12s:      %.3f     %.3f      %.3f\n',...
        cname{i},std(Tols(c,:)), std(Tswe0(c,:)),  std(Tswe1(c,:)));
end

Ta=tinv(1-alph,N-P);
fprintf('\nFPR (nominal %g, CI [%.4f,%.4f])\n',alph,alph+[-1,1]*sqrt(alph*(1-alph)/Nelm));
fprintf('                T_ols   T_swe0  T_swe0wb  T_swe1  T_swe1wb\n');
for i = 1:length(cI)
    c=cI(i);
    fprintf('%12s:   %.3f    %.3f    %.3f    %.3f    %.3f\n',...
        cname{i},mean(Tols(c,:)>=Ta),mean(Tswe0(c,:)>=Ta),mean(Pwb0(c,:)<=alph),mean(Tswe1(c,:)>=Ta),mean(Pwb1(c,:)<=alph));
end

figure;imagesc(X)

Str={'Pure between','Pure within','Mixed'};
for i = 1:length(cI)
    c=cI(i);
    figure
    P0 = tcdf(Tswe0(c,:),N-P,'upper');
    P1 = tcdf(Tswe1(c,:),N-P,'upper');
    loglog((1:Nelm)'/Nelm,[sort(P0'),sort(P1'),sort(Pwb0(c,:)'),sort(Pwb1(c,:)')],'linewidth',2);
    xlabel('Expected Quantile   (conservative above identity, invalid below)')
    ylabel('Observed Quantile')
    grid on;legend('swe0','swe1','sweWB0','sweWB1','AutoUpdate','off')
    set(refline(1),'linestyle',':')
    title(sprintf("Parameter %d: %s",c,cname{i})); % Str{i}
end


                                       
