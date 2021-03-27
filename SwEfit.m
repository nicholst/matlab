function [cbetahat, cbetaSE] = SwEfit(X,bID,Y,con)
% FUNCTION [cbetahat, cbetaSE] = SwEfit(X,bID,Y,con)
%
% Estimate betahat with OLS and obtain standard errors using the 'classic'
% Sandwich estimator.  
%
%   X   - Design matrix, N x P
%   bID - Block IDs, identifying clusters, might be family, site, subject (for repeated mes) 
%   Y   - Data, N x Nelm
%   con - Matrix of t contrasts, Ncon x P; if omitted defaults to eye(P)
%
%
% T. Nichols 24 March 2021
% See https://github.com/nicholst/matlab/blob/master/LICENSE


%
% Check arguments
%
if nargin < 3
    error('Insufficient arguments')
end
N = size(X,1);
P = size(X,2);
if ~all(size(bID)==[N,1])
    error('Block IDs bID not a N-vector')
end
Nelm = size(Y,2);
if size(Y,1)~=N
    error('Data Y not a matrix with N rows')
end
if nargin < 4
    con=eye(P);
else
    if size(con,2)~=P
        error('Contrast matrix con not a matrix with P columns')
    end
end
Ncon = size(con,1);


%
% OLS fit
%

pX      = pinv(X);
bh      = pX*Y;
res     = Y-X*bh;

%
% Compute robust (sandwich) variance estimator
%

IDs    = unique(bID)';

Bread  = pX*pX';
BreadX = Bread*X';
S      = zeros(P,P,Nelm);
E      = zeros(1,P,Nelm);

for b = IDs
    I    = (b==bID);
    Ns   = sum(I);
    % BreadX times half of Meat
    E    = BreadX(:,I)*res(I,:);
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



    