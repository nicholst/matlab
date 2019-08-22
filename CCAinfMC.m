function CCAinfMC(N,Ny,Nx,Nz,Npca,nR,nP,HuhJhun)
% Monte Carlo evaluation of CCA inference
% CCAinfMC(N,Ny,Nx,Nz,Npca,nR,nP,HuhJhun)
%
% N     - Number of subjects
% Ny    - Number of variables in Y
% Nx    - Number of variables in X
% Nz    - Number of nuisance variables (in Z)
% Npca  - Number of PCA dims to retain (empty for no PCA)
% nR    - Number of realizations
% nP    - Number of permutations per realization
% HuhJhun
%       - Use the Huh-Jhun method? (true/false)
%
% All arguments required.  Summary of results printed; nothing saved.
%
% An interecept is always used; Nz=1 means one nuisance variable (in
% addition to the intercept).
%
%----------------------------------------------------------------------------
% Anderson Winkler & Thomas Nichols
% August 2019

FreemanLane=false;

fprintf('Simulation parameters: ')
fprintf('N: %d, Ny: %d, Nx: %d, Nz: %d, Npca: %d, nR: %d, nP: %d, HJ: %d\n',...
    N,Ny,Nx,Nz,Npca,nR,nP,HuhJhun);

% For each realization:
for r = 1:nR
    fprintf('Realization %d:\n',r)
    
    % Create some random data. Use rng for repeatability:
    if isoctave
        rand('state','reset'); %#ok
    else
        rng('shuffle');
    end
    Y = randn(N,Ny);
    X = randn(N,Nx);
    Z = [randn(N,Nz-1) ones(N,1)];
    
    % Residual forming matrix:
    pZ = pinv(Z);
    Rz = eye(N) - Z*pZ;
    if HuhJhun
        % % Tom's version
        % [Q,S]  = svd(Rz);
        % S      = diag(S) < 10*eps;
        % Q(:,S) = [];
        
        % Huh Juhn's / Anderson's version
        [Q,D]  = schur(Rz);
        D      = abs(diag(D)) < 10*eps;
        Q(:,D) = [];
    else
        Q = eye(N);
    end
    
    % Remove nuisance, possibly project, make sure nothing is rank deficient:
    Y = Q'*Rz*Y;
    X = Q'*Rz*X;
    if ~ isempty(Npca)
        Y = epca(Y,Npca);
        X = epca(X,Npca);
    end
    
    % CCA:
    %[A,B,R,U,V,stats] = canoncorr(Y,X);
    [R,lW] = cca(Y,X);

    % Permutation test:
    pperm  = ones(size(R));
    pcorr  = ones(size(R));
    for p = 1:(nP-1)
        if r == 1 && p == 1
            ppermrep  = zeros(nR,length(R));  % Permutation p-value based on F-test or Wilks'
            pcorrrep  = zeros(nR,length(R));  % Permutation p-value based on r
            corrFirst = zeros(nR,length(R));  % To store the CCs for the first permutation
            corrLast  = corrFirst;            % To store the CCs for the last permutation
                                              % (any perm would do, the last is simpler as it
                                              % it stays in the memory at the end of the
                                              % permutation loop).
        end
        if FreemanLane
            tmpX = X(randperm(size(X,1)),:);
            %[~,~,Rp,~,~,statsp] = canoncorr(Y,tmpX-Z*(pZ*tmpX));
            [Rp,lWp] = cca(Y,tmpX-Z*(pZ*tmpX));
        else
            %[~,~,Rp,~,~,statsp] = canoncorr(Y,X(randperm(size(X,1)),:));
            [Rp,lWp] = cca(Y,X(randperm(size(X,1)),:));
        end
        %pperm  = pperm  + (statsp.F >= stats.F);
        pperm  = pperm  + (lWp >= lW);
        pcorr  = pcorr  + (Rp  >= R);
    end
    pperm  = pperm/nP;
    pcorr  = pcorr/nP;
    ppermrep (r,:) = pperm;
    pcorrrep (r,:) = pcorr;
    corrFirst(r,:) = R;
    corrLast (r,:) = Rp;
    fprintf('- P-values (permutation [Wilks]):'); disp(pperm);
    fprintf('- P-values (permutation [corr]):'); disp(pcorr);
    fprintf('- CCs (not permuted):');  disp(R);
    fprintf('- CCs (random perm):'); disp(Rp);
end

alpha = 0.05;
fprintf('Results:\n')
fprintf('FPR (permutation [Wilks]):');  disp(mean(ppermrep <= alpha));
fprintf('FPR (permutation [corr]):');  disp(mean(pcorrrep <= alpha));
fprintf('Mean CCs (not permuted):');  disp(mean(corrFirst));
fprintf('Mean CCs (random perm):'); disp(mean(corrLast));

% =================================================================
function [cc,lW] = cca(X,Y)
[Qx,~,~] = qr(X,0);
[Qy,~,~] = qr(Y,0);
k  = min(size(X,2),size(Y,2));
[~,D,~] = svd(Qx'*Qy,0);
cc = min(max(diag(D(:,1:k))',0),1);
lW = -fliplr(cumsum(fliplr(log(1-cc.^2))));

% =================================================================
function y = isoctave
persistent isoct;
if isempty(isoct),
    isoct = exist('OCTAVE_VERSION','builtin') ~= 0;
end
y = isoct;