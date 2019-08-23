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
    R = cca(Y,X);

    % Permutation test:
    pcorr  = ones(size(R));
    for p = 1:(nP-1)
        if r == 1 && p == 1
            pcorrrep  = zeros(nR,length(R));  % Permutation p-value based on r
            corrFirst = zeros(nR,length(R));  % To store the CCs for the first permutation
            corrLast  = corrFirst;            % To store the CCs for the last permutation
                                              % (any perm would do, the last is simpler as it
                                              % it stays in the memory at the end of the
                                              % permutation loop).
        end
        if FreemanLane
            tmpX = X(randperm(size(X,1)),:);
            Rp = cca(Y,tmpX-Z*(pZ*tmpX));
        else
            Rp = cca(Y,X(randperm(size(X,1)),:));
        end
        pcorr  = pcorr  + (Rp  >= R);
    end
    pcorr  = pcorr/nP;
    pcorrrep (r,:) = pcorr;
    corrFirst(r,:) = R;
    corrLast (r,:) = Rp;
    fprintf('- P-values (permutation [corr]):'); disp(pcorr);
    fprintf('- CCs (not permuted):');  disp(R);
    fprintf('- CCs (random perm):'); disp(Rp);
end

alpha = 0.05;
fprintf('Results:\n')
fprintf('FPR (permutation [corr]):');  disp(mean(pcorrrep <= alpha));
fprintf('Mean CCs (not permuted):');  disp(mean(corrFirst));
fprintf('Mean CCs (random perm):'); disp(mean(corrLast));

fprintf('AWK(FPR-pr,mr,mpr): %g %g %g\n',...
	[pcorrrep(:,1),...
	 corrFirst(:,1),...
	 corrLast(:,1)]');

% =================================================================
function cc = cca(X,Y)
[Qx,~,~] = qr(X,0);
[Qy,~,~] = qr(Y,0);
cc = svds(Qx'*Qy,1);

% =================================================================
function y = isoctave
persistent isoct;
if isempty(isoct),
    isoct = exist('OCTAVE_VERSION','builtin') ~= 0;
end
y = isoct;
