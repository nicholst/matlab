%
%  Demonstration of missing-variable GLM permutation, where 
%     Y is a N x K matrix of responses with irregular missing values that 
%          satisfy Missing Completely at Random,
%     X is a shared N x P matrix,
%     M is a N x K mask matrix, 1 indicating non-missingness, 0 missing
%     c is a 1 x P contrast, which in this simplfied code must have a single
%          non-zero element

Nperm = 1000;

% Invent some data
N=30;
K=500;
P=2;
Pmiss=0.5;

X = [ones(N,1) rand(N,P-1)];
Y = randn(N,K)*10;
M = rand(N,K)>Pmiss; % 'Mask' 1-present, 0-missing
c = [0 1];


% Call glm_miss once, with correctly labeled data
[z0,t0,bh0,s20] = glm_miss(X,Y,M,c);


% Permute!
zperm = zeros(Nperm,K);
Mzperm = zeros(Nperm,1);
for i = 1:Nperm
    % Use Draper-Stoneman method (Winkler et al. 2014)
    % We assume contrast c has just a single non-zero element
    Ic = c~=0;
    Xc = X(:,Ic);
    Xn = X(:,~Ic);
    P = randperm(N);
    zperm(i,:) = glm_miss([Xc(P,:) Xn],Y,M,[c(Ic) c(~Ic)]);
    Mzperm(i) = max(zperm(i,:));
end

% Compute P-values
Puncr = zeros(K,1);
Pcorr = zeros(K,1);
for k = 1:K
    Puncr(k) = (sum(z0(k)<=zperm(:,k))+1)/(Nperm+1);
    Pcorr(k) = (sum(z0(k)<=Mzperm)+1)/(Nperm+1);

end

    

