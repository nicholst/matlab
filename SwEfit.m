function [cbetahat, cbetaSE, Vg_out] = SwEfit(X,bID,Y,con,Vg)
% FUNCTION [cbetahat, cbetaSE, Vg_out] = SwEfit(X,bID,Y[,con,Vg])
%
% Estimate betahat with OLS and obtain standard errors using the 'classic'
% Sandwich estimator, 
%
%   X   - Design matrix, N x P
%   bID - Block IDs, identifying clusters, might be family, site, subject (for repeated mes) 
%   Y   - Data, N x Nelm
%   con - Matrix of t contrasts, Ncon x P; if omitted defaults to eye(P)
%   Vg  - Global working covariance; a cell array such that inv(Vg{i}) is an 
%         estimate of the covariance of the data for block i. If omitted, independence
%         working covariance is used. If is a non-zero scalar, is estimated
%
%   Vg_out  - Global averaged covariance; an estimate suitable for re-calling function in a 
%         second pass.
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
bNs    = zeros(1,Nblock);
for i = 1:Nblock
    bNs(i) = sum(IDs(i)==bID);
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
if nargin < 5
    Vg = {};
else
    if ~iscell(Vg) || length(Vg)~=Nblock
        error('Global variance not a cell array with an matrix per block');
    end
    for i = 1:Nblock
        I    = find(IDs(i)==bID);
        Ns   = length(I);
        if ~all(size(Vg{i})==[Ns Ns])
            error(sprintf('Global variance for block %d wrong size',i));
        end
        W{i} = inv(Vg{i});
    end
end
if nargout>2 
    Vg_out = cell(Nblock,1);
else
    Vg_out = cell(0);
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
        I    = find(IDs(i)==bID);
        if ~all(size(Vg{i})==[bNs(i) bNs(i)])
            error(sprintf('Global variance for block %d wrong size',i));
        end
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
    I     = find(IDs(i)==bID);
    Ra{i} = sqrtm(inv(eye(bNs(i))-X(I,:)*BreadXW(:,I)));
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
    I    = find(IDs(i)==bID);
    % BreadXW times half of Meat, with type 1 residual correction
    ares = Ra{i}*res(I,:);
    E    = BreadXW(:,I)*ares;
    E    = reshape(E,[P,1,Nelm]);
    % Full `Bread*Meat*Bread' contribution for block b
    S    = S + pagemtimes(E,'none',E,'transpose');
    if length(Vg_out)>0
        Vg_out{i} = ares*ares'/Nelm;
    end
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



    