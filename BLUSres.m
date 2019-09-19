function [u,I1,I0] = BLUSres(Y,X,varargin)
% Compute BLUS (Best Linear Unbiased Scalar covariance) residuals
% FORMAT [u,I1,I0] = BLUSres(Y,X[,I1])
%
% Y  - Data  - May have multiple columns
% X  - Design - Must have same number of rows as Y
% I1 - Subset of indicies retained; if omitted randomly computed as
%      I=randperm(N,N-P)
%
% u  - BLUS residuals
% I1 - Indicies retained
% I0 - Indicies dropped
%
%_______________________________________________________________________
% TE Nichols Sept 2019

N = size(Y,1);
P = size(X,2);
if size(X,1)~=N
  error('Rows and and X and Y don''t match')
end

if (nargin >= 3)
  I1 = varargin{1};
  I0 = setdiff(1:N,I1);
  if rank(X(I0,:))<P
    error('Omitted subset not full rank')
  end
else
  Done=false
  while (~Done)
    I1 = randperm(N,N-P);
    I0 = setdiff(1:N,I1);
    if rank(X(I0,:))==P
      Done=true;
    end
  end
end

% Magnus, J. R., & Sinha, A. K. (2005). On Theil’s errors. The Econometrics
% Journal, 8(1), 39–54.

M    = eye(N) - X*inv(X'*X)*X';
S    = eye(N); 
S    = S(:,I1);
A    = M*S*sqrtm(inv(S'*M*S));

u    = A'*Y;

% You can verify that A'*X = 0 and A'*A = I_{N-P}

return

%
% Multiple versions that work or don't
%

% Vinod, H. D. (2014). Theil’s BLUS Residuals and R Tools for Testing and 
% Removing Autocorrelation and Heteroscedasticity. SSRN Electronic Journal,
% (1), 1–9. https://doi.org/10.2139/ssrn.2412740 

% Works, but more complicated; agrees with Magnus & Sinha

X0 = X(I0,:);
X1 = X(I1,:);

r  = Y - X*pinv(X)*Y;

r0 = r(I0,:);
r1 = r(I1,:);

[V,D] = eig(X0*inv(X'*X)*X0');
D     = diag(D);
Z     = zeros(P,P);
for i = 1:P
  Z = Z + sqrt(D(i))/(1+sqrt(D(i)))*V(:,i)*V(:,i)';
end

u     = r1 - X1*inv(X0)*Z*r0;


% Grenier, M., Léger, C., Grenier, M., & Leger, C. (2000). Bootstrapping
% Regression Models with BLUS Residuals. The Canadian Journal of 
% Statistics, 28(1), 31. https://doi.org/10.2307/3315880

% Completely different from Magnus and Sinha

M     = eye(N) - X*inv(X'*X)*X';
M11   = (M(I1,I1)+M(I1,I1)')/2;
[V,D] = eig(M11);
A     = [ -V*sqrt(D)*V'*X1*X0 V*sqrt(D)*V' ];

u     = A*Y;


% Linear Models and Time-Series Analysis: Regression, ANOVA, ARMA and GARCH
% By Marc S. Paolell, Section 1.5, pp 47-50.

% Close but not quite same as Magnus & Sinha.

Za   = X1*inv(X0); 
D    = eig(X0*inv(X'*X)*X0');
i1   = find(D<1 & D>0); 
ni   = size(i1,1);
D    = [D(i1);ones(N-P-ni,1)]; 
D    = sort(D); 
D    = diag(D);
[V, tempD] = eig(eye(N-P) + Za*Za');
tempD= diag(tempD); 
[tempD, i2] = sortrows(tempD);
V    = V(:,i2(end:-1:1)); 
C_1  = V*D*V'; 
C    = [-C_1*Za C_1];
u2    = C*Y;


