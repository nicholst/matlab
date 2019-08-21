CCAinfMC
% Monte Carlo evaluation of CCA inference
% CCAinfMC
%
% Following varibles must be set in the workspace in order for this to run
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
% An interecept is always used; Nz=1 means one nuisance variable (in
% addition to the intercept).
%
%----------------------------------------------------------------------------
% Anderson Winkler & Thomas Nichols
% August 2019

FreemanLane=false

% Vars for later:
ppararep = zeros(nR,1);  % Parametric p-value from CCA F-test
ppermrep = zeros(nR,1);  % Permutation p-value based on F-test
pcorrrep = zeros(nR,1);  % Permutation p-value based on r
corrrep  = zeros(nR,1);  % First canonical correlation
corrrep0m= zeros(nR,1);  % Mean of permuted first canonical correlation

% For each realization:
for r = 1:nR
    fprintf('%d\n',r)
    
    % Create some random data. Use rng for repeatability:
    rng('shuffle');
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
    Y   = Q'*Rz*Y;
    if ~isempty(Npca)
      Y   = epca(Y,Npca);
    end
    X   = Q'*Rz*X;
    if ~isempty(Npca)
      X   = epca(X,Npca);
    end
    
    % CCA:
    [A,B,R,U,V,stats] = canoncorr(Y,X);
    
    % Permutation test:
    pperm = ones(size(R));
    pcorr = ones(size(R));
    corrs = zeros(nP-1,1);
    for p = 1:(nP-1)
        if FreemanLane
	  tmpX=X(randperm(size(X,1)),:);
	  [~,~,Rp,~,~,statsp] = canoncorr(Y,tmpX-Z*(pZ*tmpX));
	else
	  [~,~,Rp,~,~,statsp] = canoncorr(Y,X(randperm(size(X,1)),:));
	end
        pperm = pperm + (statsp.pF <= stats.pF);
        pcorr = pcorr + (Rp >= R);
	corrs(p) = Rp(1);
    end
    pperm = pperm/nP;
    pcorr = pcorr/nP;
    ppararep(r) = statsp.pF(1);
    ppermrep(r) = pperm(1);
    pcorrrep(r) = pcorr(1);
    corrrep(r)  = R(1);
    corrrep0m(r)= mean(corrs);
end

fprintf('P,nPF,nPr,r,r0m: %g %g %g %g %g\n', ...
	[ppararep ppermrep pcorrrep corrrep corrrep0m]')


