function [cbetahat, cbetaSE, Vg] = SwEfit(X,bID,Y,con,Vg,RA)
% FUNCTION [cbetahat, cbetaSE, Vg] = SwEfit(X,bID,Y[,con,Vg,RA])
%
% Estimate betahat with OLS and obtain standard errors using the 'classic'
% Sandwich estimator, 
%
%   X   - Design matrix, N x P
%   bID - Block IDs, identifying clusters, might be family, subject (for repeated meas) 
%   Y   - Data, N x Nelm
%   con - Matrix of t contrasts, Ncon x P; if omitted or empty defaults to eye(P)
%   Vg  - Global working covariance:
%             []  - use independence (default)
%             1   - estimate the global covariance on the fly, otw
%             Vg  - cell array such that Vg{i} is the working covariance for block i
%
%   Vg  - Global working covariance; useful if Vg is computed on the fly.
%   RA  - Residual adjustment; options are
%            'HC0' - No adjustment
%            'HC1' - Scalar adjustment, sqrt(N/(N-P))
%            'HC2' - Diagonal adjustment, 1/sqrt(1-hii)
%             'C2' - Multivariate HC2, (1-Hii)^(-1/2), where H_ii is the sub-matrix of
%                    the hat matrix H.
%            'HC3' - Diagonal adjustment, extra correction, 1/(1-hii)
%            'HC4' - Diagonal adjustment, extra EXTRA correction, 1/(1-hii)^(delta/2),
%                    where delta = min(4,hii/mean(hii))
%
% T. Nichols 11 Jan 2022
% See https://github.com/nicholst/matlab/blob/master/LICENSE


%
% Check arguments & set basic variables
%
if nargin < 3
    error('Insufficient arguments')
end
N = size(X,1);
P = size(X,2);
if ~all(size(bID)==[N,1])
    error('Block IDs bID not a N-vector')
end
IDs    = unique(bID(:));
Nblock = length(IDs);
bI     = cell(Nblock,1);
bN     = zeros(1,Nblock);
for i = 1:Nblock
    bI{i} = find(IDs(i)==bID);
    bN(i) = length(bI{i});
end
Nelm   = size(Y,2);
if size(Y,1)~=N
    error('Data Y not a matrix with N rows')
end
if nargin < 4 || isempty(con)
    con=eye(P);
else
    if size(con,2)~=P
        error('Contrast matrix con not a matrix with P columns')
    end
end
Ncon = size(con,1);
CalcVg = false;
if nargin < 5
    Vg = [];
else
    if ~iscell(Vg) 
        if ~isempty(Vg) && prod(size(Vg))~=1
            error('Global variance Vg not a scalar nor a cell array with an matrix per block')
        end
        if Vg
            CalcVg = true;
            Vg = cell(Nblock,1);
        end
    else
        if length(Vg)~=Nblock
            error('Global variance Vg not a scalar nor a cell array with an matrix per block');
        end
        for i = 1:Nblock
            if ~all(size(Vg{i})==[bN(i) bN(i)])
                error(sprintf('Global variance for block %d wrong size',i));
            end
        end
    end
end
if nargin < 6
    RA = 'C2';
else
    if ~any(strcmp(lower(RA),lower({'HC0','HC1','HC2','C2','HC3','HC4','HC4m','HC5'})))
        error('Unknown residual adjustment option')
    end
end

%
% Compute global working covariance, if needed
%
if CalcVg
    pX  = pinv(X);
    res = Y-X*pX*Y;
    for i = 1:Nblock
        I    = bI{i};
        % 'C2' adjustment, a `block-wise 1/sqrt(1-hii)'
        Ra   = sqrtm(inv(eye(bN(i))-X(bI{i},:)*pX(:,I)));
        ares = Ra*res(I,:);
        Vg{i} = ares*ares'/Nelm;
    end
end

%
% Compute BreadXW = inv(X'*W*X)*X'*W, used to estimate bh and S
%

if length(Vg)==0
    BreadXW = pinv(X);
else
    XtWX = zeros(P,P);
    XtW  = zeros(P,N);
    XtWh = zeros(P,N);
    for i = 1:Nblock
        I        = bI{i};
        W{i}     = pinv(Vg{i});
        Wh{i}    = sqrtm(W{i});
        XtWX     = XtWX + X(I,:)'*W{i}*X(I,:);
        XtW(:,I) = X(I,:)'*W{i};
        XtWh(:,I)= X(I,:)'*Wh{i};
    end
    BreadXW  = pinv(XtWX)*XtW;
    BreadXWh = pinv(XtWX)*XtWh;
end

%
% Residual adjustment
%
% Usual robust standard error name conventions; see
% https://cran.r-project.org/package=sandwich/vignettes/sandwich.pdf
%
% Note, when using a (non identity) working covariance hii values can be
% less than 0 (normally, 0 <= hii <= 1).
%
Ra = cell(1,Nblock);
H  = XtWh'*BreadXWh;
mH = mean(diag(H));
mxH= max(diag(H));
for i = 1:Nblock
    I     = bI{i};
    Hii   = H(I,I);
    hii   = diag(H(I,I));
    switch RA
      case 'HC0'
        Ra{i} = 1;
      case 'HC1'
        Ra{i} = sqrt(N/(N-P));
      case 'HC2'
        Ra{i} = diag((1-hii).^(-0.5));
      case 'C2'
        Ra{i} = sqrtm(inv(eye(bN(i))-Hii));
      case 'HC3'
        Ra{i} = diag((1-hii).^(-1));
      case 'HC4'
        % TN modification: Original HC4 uses a delta of
        %  min(4,hii/mH))
        % here I ensure delta never falls below mH/10
        delta = max(min(4,hii/mH),mH/10);
        Ra{i} = diag((1-hii).^(-delta/2));
      case 'HC4m'
        delta = min(1,hii/mH)+min(1.5,hii/mH);
        Ra{i} = diag((1-hii).^(-delta/2));
      case 'HC5'
        % TN modification: Original HC5 uses a delta of
        %   min(max(4,0.7*mxH/mH),hii/mH)
        % but on horribly skewed data that causes crazy large powers
        % Here, I put a cap of 10, which seems reasonable, and
        % bound below by mH/10
        delta = max(min(max(4,min(10,0.7*mxH/mH)),hii/mH),mH/10);
        Ra{i} = diag((1-hii).^(-delta/2));
    end
end

%
% OLS (or working-covariance modified) fit
%

bh   = BreadXW*Y;
res  = Y-X*bh;


%
% Compute robust (sandwich) variance estimator
%
S      = zeros(P,P,Nelm);
for i = 1:Nblock
    I    = bI{i};
    % BreadXW times half of Meat, with type 1 residual correction
    ares = Ra{i}*res(I,:);
    E    = BreadXW(:,I)*ares;
    E    = reshape(E,[P,1,Nelm]);
    % Full `Bread*Meat*Bread' contribution for block b
    S    = S + pagemtimes(E,'none',E,'transpose');
end

%
% Compute contrasts
%
cbetahat = zeros(Ncon,Nelm);
cbetaSE  = zeros(Ncon,Nelm);

for i = 1:Ncon
    c             = con(i,:);
    cbetahat(i,:) = c*bh;
    cbetaSE(i,:)  = sqrt(pagemtimes(pagemtimes(c,S),'none',c,'transpose'));
end


end

