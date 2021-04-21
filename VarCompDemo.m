N = 100;     % Number of subjects
M = 200;    % Number of voxels/imaging features combined in generative modoel
Mpca = 20;  % For linear regression approach, number of pca components to use
nRlz = 1000; % Number of releasations

sig2M = 0.4; % Morphometric variance -- Set this [0,1] so that m2 is same as sig2M
sig2E = 1-sig2M;
m2 = sig2M / (sig2M + sig2E); % Morphometricity


% Morphometry 
X = ones(N,1); beta = [1];
Z = rand(N,M) + 0.1*linspace(0,1,N)';
Z = Z - mean(Z); % Center each voxel over subjects
Z = Z ./ std(Z); % Standardize each voxel over subjects

% Mophmetric similarity -- Compute similarity between subjects over voxels
% (This is *not* the correlation between subjects -- that would require centering
% standardising each subject's data, instead of each voxel).
B = Z*Z'/M; 

[V,D]=eig(B);
D1 = [diag(D) ones(N,1)];

% XB - PCA based approximation of Z
if Mpca>=M
    Mpca=M;
    XB = Z;
else
    XB = V(:,end-Mpca+1:end);
end

% Plot
figure(1)
imagesc(B);title('NxN intersubject similarity')
figure(2)
plot(diag(D));title('Eigenvalues of B')

% Simulate
[sig2Mh,sig2Eh,R2XB]=deal(zeros(nRlz,1)); % sig2M-hat, sig2E-hat
for i=1:nRlz
    u   = randn(M,1)*sqrt(sig2M/M);
    eps = randn(N,1)*sqrt(sig2E);

    Y = X*beta + Z*u + eps;
    
    % Estimate variance components with FPHI approach
    % Ganjgahi et al. (2015) https://doi.org/10.1016/j.neuroimage.2015.03.005 
    Ys = V'*Y;
    Xs = V'*X;
    es = Ys - Xs*pinv(Xs)*Ys;
    tmp = pinv(D1)*(es.^2);
    sig2Mh(i) = tmp(1);
    sig2Eh(i) = tmp(2);

    % Estimate linear regression approximation
    bbh = pinv([X XB])*Y;  % "beta_b - hat"
    SStot = sum((Y - mean(Y)).^2);
    SSexp = sum(([XB]*bbh(size(X,2)+1:end)).^2);
    R2XB(i) = SSexp/SStot;
end

% Estimate m2-hat
m2h = sig2Mh ./ (sig2Mh + sig2Eh);


% Report...

fprintf('N: %d   M: %d   Mpca: %d   nRlz: %d   \n',N,M,Mpca,nRlz);
fprintf('Tru-m2: %0.4f  Est-m2: %0.4f  R2_XB: %0.4f\n',...
        m2,mean(m2h),mean(R2XB))

% Gory detail on bias/variance/mse... needed to make sure I was estimating 
% variaince components well
if (0)
    RelBias=@(x,mu)[(mean(x)-mu)/mu];
    RelStd =@(x,mu)[std(x)/mu];
    RelrMSE=@(x,mu)[sqrt(mean((x-mu).^2))/mu]; % relative root MSE
    fprintf('sig2M: Bias % .1f%%  Std %5.1f%%  rMSE %5.1f%%\n',...
            100*RelBias(sig2Mh,sig2M),...
            100*RelStd( sig2Mh,sig2M),...
            100*RelrMSE(sig2Mh,sig2M));
    fprintf('sig2E: Bias % .1f%%  Std %5.1f%%  rMSE %5.1f%%\n',...
            100*RelBias(sig2Eh,sig2E),...
            100*RelStd( sig2Eh,sig2E),...
            100*RelrMSE(sig2Eh,sig2E));
    fprintf('m2:    Bias % .1f%%  Std %5.1f%%  rMSE %5.1f%%\n',...
            100*RelBias(m2h,m2),...
            100*RelStd( m2h,m2),...
            100*RelrMSE(m2h,m2));
end
