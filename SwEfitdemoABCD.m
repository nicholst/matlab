% SwEdemo code
% Attempting to find minimum run time for mass-univariate application of sandwich estimator
%
% Based on https://github.com/nicholst/tenR/blob/master/SwEdemo.R
%
% T. Nichols 13 Jan 2022
% See https://github.com/nicholst/matlab/blob/master/LICENSE

% 'block' is cluster variable, might be family, site, subject (for repeated mes)
Nelm      = 1000; % Number of vertices/voxels
rho       = 0.95; % Intrablock correlation... maxed out to verify SwE is working
alph      = 0.05; % Nominal alpha for FPR report
nWB       = 100;    % Number of Wild Bootstrap iterations - Start at 100; try 1000

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

% AMD: changed to handle design matrix in original form  
try
  X = load('FamSubVisX.dat');
  Nms = readcell('FamSubVisX_names.txt', 'Delimiter',"");
catch
  tbl = readtable('ABCD_rel3.0_long_desmat_PCs_SES_interview_age.txt');
  reordervec = [3 1 2];
  X = NaN(size(tbl));
  X(:,4:end) = table2array(tbl(:,4:end));
  for j = 1:length(reordervec)
    tmp = tbl(:,reordervec(j));
    [dummy IA IC] = unique(tmp,'stable');
    X(:,j) = IC;
  end
  Nms = tbl.Properties.VariableNames;
end

Fam   = X(:,1);
Sub   = X(:,2);
Vis   = X(:,3);
Age   = X(:,4);
X     = X(:,5:end);
Nms(1:4)=[];

Iblock = Fam;
Nperblock = groupcounts(Iblock);

% Need rows in family order to simplify simulation (also sort so largest 
% families are at the end, easier to see block structure)
[~,I]  = sort(Iblock);
[~,II] = sortrows([Iblock(I),repelem(Nperblock,Nperblock)],[2 1]);
X      = X(I(II),:);
Fam    = Fam(I(II),:);
Sub    = Sub(I(II),:);
Vis    = Vis(I(II),:);
Age    = Age(I(II),:);
Iblock = Iblock(I(II));
Nperblock=sort(Nperblock);

Opt=''
%Vg={    0,   0,    0,    0,     0,    0,    0,     0,    1,   1,    1,    1,    1,    1,     1,    1}; % Global estimated working cov
%RA={'HC2','C2','HC3','HC4','HC4m','HC5','HC5m','HC6','HC2','C2','HC3','HC4','HC4','HC5','HC5m','HC6'}; % Different res adj methods
Vg={    0,    0,    1,    1,   1}; % Global estimated working cov
RA={'HC5','HC6','HC3','HC5','HC6'}; % Different res adj methods
nSwE = length(Vg);
cI=[1,2,6,10,15,47,53,62];
switch Opt
  case 'PureWithin'
    for s = unique(Sub)'
        If = Sub==s;
        X(If,2) = X(If,2) - mean(X(If,2));
    end
  case 'GaussCov'
    X(:,3:end-1) = randn(N,size(X,2)-3);
    Nms(3:end-1) = {'rand'};
  case 'BalBinCov'
    X(:,3:end-1) = binor(1,0.5,N,size(X,2)-3);
    Nms(3:end-1) = {'rand'};
  case 'NoCov'
    X=X(:,[1:2,end]);
    Nms=Nms([1:2,end]);
  case 'BadCovOnly'
    I=[1,2,20,53,size(X,2)];
    X=X(:,I);
    Nms=Nms(I);
    cI=[2,3,4];
  case 'FilterByFamilySize'
    MinSz=2
    I      = repelem(Nperblock,Nperblock)>=MinSz;
    X      = X(I,:);
    Fam    = Fam(I,:);
    Sub    = Sub(I,:);
    Vis    = Vis(I,:);
    Age    = Age(I,:);
    Iblock = Iblock(I);
    Nperblock=Nperblock(Nperblock>=MinSz);
    disp(sprintf('Only using familes of %d or more (%d rows)',MinSz,sum(I)))
end

N      = size(X,1);
P      = size(X,2);
Nblock = length(unique(Fam));

% AMD: generate synthesized data with F and S random effects
if 0
  % Simulate repeated measures data, N x Nelm, with intrablock correlation rho
  Y = sqrt(1-rho)*randn(N,Nelm) + ...
      sqrt(rho)  *repelem(randn(Nblock,Nelm),Nperblock,1);
else
  sig2vec = [0.5 0.5]*0.95; % Half F and S
  Y = sqrt(1-sum(sig2vec))*randn(N,Nelm);
  for fi = 1:max(Fam)
    ivec = find(Fam==fi);
    Y(ivec,:) = Y(ivec,:) + sqrt(sig2vec(1))*randn(1,Nelm);
  end
  for si = 1:max(Sub)
    ivec = find(Sub==si);
    Y(ivec,:) = Y(ivec,:) + sqrt(sig2vec(2))*randn(1,Nelm);
  end
end

% OLS fit
tic;
pX      = pinv(X);
hii     = diag(X*pX);
bh      = pX*Y;
res     = Y-X*bh;
sig2ols = sum(res.^2)/(N-P);
SEols   = sqrt(diag(pX*pX') .* sig2ols);
Tols    = bh./SEols;
fprintf('%-20s','OLS');toc
[SwElab,cbetahat,cbetaSE,Tswe,Pswe] = deal(cell(0));

jj=1;
SwElab{jj}='OLS';
cbetahat{jj}=bh;
cbetaSE{jj}=SEols;
Tswe{jj}=Tols;
Pswe{jj}= tcdf(Tols, N-P,'upper');

for j=1:nSwE
    jj=jj+1;
    SwElab{jj} = sprintf('SwE:Vg%d,%s',Vg{j},RA{j});
    tic;
    [cbetahat0,cbetaSE0]=SwEfit(X,Iblock,Y,[],Vg{j},RA{j});
    fprintf('%-20s',SwElab{jj});toc
    Tswe0=cbetahat0./cbetaSE0;
    cbetahat{jj}=cbetahat0;
    cbetaSE{jj}=cbetaSE0;
    Tswe{jj}=Tswe0;
    Pswe{jj}=tcdf(Tswe0,N-P,'upper');

    if nWB>0
        jj=jj+1;
        SwElab{jj} = sprintf('SwE:Vg%d,%s,WB',Vg{j},RA{j});
        tic
        % Wild Bootstrap
        res0 = (Y-X*cbetahat0)./sqrt(1-hii);
        % Cribari-Neto et al. (2004) suggests that a HC4-style standardisation 
        % works better
        Pwb0 = zeros(P,Nelm); 
        for i=1:nWB
            if rem(i,10)==0; fprintf('%d ',i); end
            WBf=repelem(2*binornd(1,0.5,Nblock,1)-1,Nperblock,1);
            Ywb0 = WBf.*res0;
            [cbwb,cbSEwb] = SwEfit(X,Iblock,Ywb0,[],Vg{j},RA{j});
            Pwb0 = Pwb0 + (cbwb./cbSEwb >= Tswe0);
        end
        Pwb0 = (Pwb0+1)/(nWB+1);
        fprintf('%20s\n ',SwElab{j});toc
        cbetahat{jj}=[];
        cbetaSE{jj}=[];
        Tswe{jj}=[];
        Pswe{jj} = Pwb0;
    end    
end
nSwE=length(Pswe);

% Compute summaries
[relSD,Tstd,FPR,FPR01]=deal(zeros(P,nSwE));
for i = 1:P
    for jj=1:nSwE
        if isempty(cbetahat{jj})
            relSD(i,jj)=NaN;
            Tstd(i,jj)=NaN;
        else
            relSD(i,jj)=std(cbetahat{1}(i,:))/std(cbetahat{jj}(i,:));
            Tstd(i,jj)=std(Tswe{jj}(i,:));
        end
        if isnan(Pswe{jj}(i,1))
            FPR(i,jj)=NaN
            FPR01(i,jj)=NaN
        else
            FPR(i,jj)=mean(Pswe{jj}(i,:)<=alph);
            FPR01(i,jj)=mean(Pswe{jj}(i,:)<=0.01);
        end
    end
end

% Select some variables
cname=Nms(cI);
for i = 1:length(cname)
    if length(cname{i})>20
        cname{i}=cname{i}(end-19:end);
    end
end

% Print tabular results
SwElab=strrep(SwElab,'SwE:','');

fprintf('\nEfficiency: SD(OLS)/SD(SwE) \n%20s ',' ')
fprintf('%s',sprintf('%12s',SwElab{:}))
for i = 1:length(cI)
    fprintf('\n%20s:%s',cname{i},sprintf('%12.3f',relSD(cI(i),:)))
end

fprintf('\n\nStandard deviation of T (should be 1.0)\n%20s ',' ');
fprintf('%s',sprintf('%12s',SwElab{:}));
for i = 1:length(cI)
    fprintf('\n%20s:%s',cname{i},sprintf('%12.3f',Tstd(cI(i),:)))
end

fprintf('\n\nFPR (nominal %g, CI [%.3f,%.3f])\n%20s ',alph,alph+[-1,1]*sqrt(alph*(1-alph)/Nelm),' ');
fprintf('%s',sprintf('%12s',SwElab{:}));
for i = 1:length(cI)
    c=cI(i);
    fprintf('\n%20s:%s',cname{i},sprintf('%12.3f',FPR(cI(i),:)))
end
fprintf('\n');

%figure;imagesc(X)

%         1 2 3 4 5 6 7 8 9 10 11
% No WB
%SwEomit=[     3   5   7   9    11];
% Only WB
SwEomit=[ 1  2  4   6   8   10   ];


for i = 1:length(cI)
    c=cI(i);
    I=[1:Nelm]';
    figure
    Ps = zeros(Nelm,0);
    for jj = 1:nSwE
        if any(SwEomit==jj)
            Ps(:,jj) = NaN;
        else
            Ps(:,jj) = sort(Pswe{jj}(c,:));
        end
    end
    h=loglog(I/Nelm,Ps);
    SetDefLines(h,2)
    xlabel('Expected Quantile   (conservative above identity, invalid below)')
    ylabel('Observed Quantile')
    grid on
    legend(SwElab{:},'AutoUpdate','off')
    set(refline(1),'linestyle',':')
    bI=linspace(1,Nelm,25)';
    rI=flipud(bI);
    xx=[bI;rI]/(Nelm+1);
    % https://en.wikipedia.org/wiki/Order_statistic#Order_statistics_sampled_from_a_uniform_distribution
    yy=[betainv(alph/2,bI,Nelm-bI+1);betainv(1-alph/2,rI,Nelm-rI+1)];
    h=get(gca,'children');
    hold on;
    fill(-log10(xx),-log10(yy),[1 1 1]*0.85,'linestyle','none','facealpha',.5)
    hold off
    uistack(h,'top');
    title(sprintf("Parameter %d: %s",c,cname{i}),'interpreter','none');
end

if (0)
figure
cols=floor(sqrt(P+1)*1.2);
rows=ceil(P/cols);
for c = 1:P
    I=[1:Nelm]';
    subplot(rows,cols,c);
    Ps = zeros(Nelm,0);
    for jj = 1:nSwE
        if any(SwEomit==jj)
            Ps(:,jj) = NaN;
        else
            Ps(:,jj) = sort(Pswe{jj}(c,:));
        end
    end
    h=loglog(I/Nelm,Ps);
    xlim(10.^[-2.5 0]);ylim(10.^[-2.5 0])
    SetDefLines(h,2);
    grid on; set(gca, 'XMinorGrid','off', 'YMinorGrid','off')
    set(refline(1),'linestyle',':')
    bI=linspace(1,Nelm,25)';
    rI=flipud(bI);
    xx=[bI;rI]/(Nelm+1);
    % https://en.wikipedia.org/wiki/Order_statistic#Order_statistics_sampled_from_a_uniform_distribution
    yy=[betainv(alph/2,bI,Nelm-bI+1);betainv(1-alph/2,rI,Nelm-rI+1)];
    h=get(gca,'children');
    hold on;
    fill(-log10(xx),-log10(yy),[1 1 1]*0.85,'linestyle','none','facealpha',.5)
    hold off
    uistack(h,'top');
    title(sprintf("%d: %s",c,...
                  extractAfter(Nms{c},max(0,length(Nms{c})-14))),...
                 'interpreter','none');
end
subplot(rows,cols,P+1)
h=loglog(xx,repmat(xx,1,nSwE));
SetDefLines(h,2);
xlabel('Expected')
ylabel('Observed')
legend(SwElab{:})
end


[aVarSkew,VarSkew,aVarKurt,VarKurt]=deal(repmat(NaN,1,P));
X0=X(1:end-1);
for i = 1:P-1
    fprintf('%d ',i)
    Xi = X(:,i);
    X0 = X(:,setdiff(1:P,i));
    Xia= Xi-X0*pinv(X0)*Xi;
    VarSkew(i) = skewness(Xi);
    VarKurt(i) = kurtosis(Xi);
    aVarSkew(i) = skewness(Xia);
    aVarKurt(i) = kurtosis(Xia);
end

figure
cols=floor(sqrt(nSwE+1)*1.2);
rows=ceil(nSwE/cols);
for jj=1:nSwE
    subplot(cols,rows,jj)
    semilogx(aVarKurt,FPR01(:,jj)/0.01,'o');abline('h',1)
    title(sprintf('RelFPR vs Kurtosis - %s',SwElab{jj}))
end
figure
for jj=1:nSwE
    subplot(cols,rows,jj)
    plot(abs(aVarSkew).^(1/2).*sign(VarSkew),FPR01(:,jj)/0.01,'o');abline('h',1)
    title(sprintf('RelFPR vs Skew - %s',SwElab{jj}))
end
figure
cols=floor(sqrt(nSwE+1)*1.2);
rows=ceil(nSwE/cols);
for jj=1:nSwE
    subplot(cols,rows,jj)
    semilogx(aVarKurt,Tstd(:,jj),'o');abline('h',1)
    title(sprintf('RelFPR vs Kurtosis - %s',SwElab{jj}))
end
figure
for jj=1:nSwE
    subplot(cols,rows,jj)
    plot(abs(aVarSkew).^(1/2).*sign(VarSkew),Tstd(:,jj),'o');abline('h',1)
    title(sprintf('RelFPR vs Skew - %s',SwElab{jj}))
end

