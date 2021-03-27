function [cbetahat, cbetaSE, Vg] = SwEfit(X,bID,Y,con,Vg)
% FUNCTION [cbetahat, cbetaSE, Vg] = SwEfit(X,bID,Y[,con,Vg])
%
% Estimate betahat with OLS and obtain standard errors using the 'classic'
% Sandwich estimator, 
%
%   X   - Design matrix, N x P
%   bID - Block IDs, identifying clusters, might be family, site, subject (for repeated mes) 
%   Y   - Data, N x Nelm
%   con - Matrix of t contrasts, Ncon x P; if omitted or empty defaults to eye(P)
%   Vg  - Global working covariance:
%             []  - use independence (default)
%             1   - estimate the global covariance on the fly, otw
%             Vg  - cell array such that inv(Vg{i}) is the working covariance for block i.
%
%   Vg  - Global working covariance; useful if Vg is computed on the fly.
%
% T. Nichols 24 March 2021
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
        end
    else
        if length(Vg)~=Nblock
            error('Global variance Vg not a scalar nor a cell array with an matrix per block');
        end
        for i = 1:Nblock
            if ~all(size(Vg{i})==[bN(i) bN(i)])
                error(sprintf('Global variance for block %d wrong size',i));
            end
            W{i} = inv(Vg{i});
        end
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
    for i = 1:Nblock
        I        = bI{i};
        W{i}     = inv(Vg{i});
        XtWX     = XtWX + X(I,:)'*W{i}*X(I,:);
        XtW(:,I) = X(I,:)'*W{i};
    end
    BreadXW = inv(XtWX)*XtW;
end

%
% Residual adjustment
%
% For block i, (I-H_ii)^-0.5, where H_ii is the sub-matrix of the hat matrix H.
% More sophisticated alternative to sqrt(n/(n-p)).
%
Ra = cell(1,Nblock);
for i = 1:Nblock
    I     = bI{i};
    Ra{i} = sqrtm(inv(eye(bN(i))-X(I,:)*BreadXW(:,I)));
end

%
% OLS fit
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



    