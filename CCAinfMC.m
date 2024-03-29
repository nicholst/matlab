function CCAinfMC(N,Ny,Nx,Nz,Npca,nR,nP,HuhJhun,sig)
% Monte Carlo evaluation of CCA inference
% CCAinfMC(N,Ny,Nx,Nz,Npca,nR,nP,HuhJhun,sig)
%
% N     - Number of subjects
% Ny    - Number of variables in Y
% Nx    - Number of variables in X
% Nz    - Number of nuisance variables (in Z)
% Npca  - Number of PCA dims to retain (empty for no PCA)
% nR    - Number of realizations
% nP    - Number of permutations per realization
% HuhJhun
%       - Use the Huh-Jhun method?
%            0 - no
%            1 - Shur
%            2 - SVD
%            3 - BLUS residuals
% sig   - Standard deviation of signal to add to X & Y (default: 0)
%
%
% If Nz<0, then Z will be derived as an average of columns of X & Y, which 
% could be one of the most pathological nuisance settings.
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

ZdepXY = false;
if (Nz<0)
  ZdepXY = true;
  Nz = abs(Nz);
end


% For each realization:
for r = 1:nR
    fprintf('Realization %d\n',r)
    
    % Create some random data. Use rng for repeatability:
    rng('shuffle');
    Y = randn(N,Ny);
    X = randn(N,Nx);
    Z = [randn(N,Nz) ones(N,1)];

    if ZdepXY
      % Add nuisance
      d = min([Nz,Nx,Ny]);
      Y(:,1:d) = Y(:,1:d) + Z(:,1:d);
      X(:,1:d) = X(:,1:d) + Z(:,1:d);
    end
    
    if sig>0
      % Add common signal
      CommSig=randn(N,1)*sig;
      Y = Y + CommSig;
      X = X + CommSig;
    end

    % Residual forming matrix:
    pZ = pinv(Z);
    Rz = eye(N) - Z*pZ;
    switch HuhJhun
     case 0
      % None
      H = Rz;
     case 1
      % Tom's version (SVD)
      [Q,S]  = svd(Rz);
      S      = diag(S) < sqrt(eps);
      Q(:,S) = [];
      H      = Q'*Rz;
     case 2
      % Huh Juhn's / Anderson's version (Schur)
      [Q,D]  = schur(Rz);
      D      = abs(diag(D)) < 10*eps;
      Q(:,D) = [];
      H      = Q'*Rz;
     case 3
      % BLUS residuals
      H      = BLUSmtx(Z);
    end
    
    % Remove nuisance, possibly project, make sure nothing is rank deficient:
    Y = H*Y;
    X = H*X;
    if ~isempty(Npca) && Npca>0
        Y = epca(Y,Npca);
        X = epca(X,Npca);
    end
    
    % CCA:
    [A,B,R,U,V,stats] = canoncorr(Y,X);

    % Permutation test:
    pperm = ones(size(R));
    pcorr = ones(size(R));
    for p = 1:(nP-1)
        if FreemanLane
            tmpX = X(randperm(size(X,1)),:);
            [~,~,Rp,~,~,statsp] = canoncorr(Y,tmpX-Z*(pZ*tmpX));
        else
            [~,~,Rp,~,~,statsp] = canoncorr(Y,X(randperm(size(X,1)),:));
        end
        pperm = pperm + (statsp.F >= stats.F);
        pcorr = pcorr + (Rp >= R);
    end
    if (r==1)
      ppararep  = zeros(nR,length(R));  % Parametric p-value from CCA F-test
      ppermrep  = zeros(nR,length(R));  % Permutation p-value based on F-test
      pcorrrep  = zeros(nR,length(R));  % Permutation p-value based on r
      corrFirst = zeros(nR,length(R));  % To store the CCs for the first permutation
      corrLast  = corrFirst;            % To store the CCs for the last permutation
					% (any perm would do, the last is simpler as it
					% it stays in the memory at the end of the
					% permutation loop).
    end
    pperm = pperm/nP;
    pcorr = pcorr/nP;
    ppararep (r,:) = stats.pF;
    ppermrep (r,:) = pperm;
    pcorrrep (r,:) = pcorr;
    corrFirst(r,:) = R;
    corrLast (r,:) = Rp;
end

fprintf('Simulation parameters:\n')
fprintf('N: %d, Ny: %d, Nx: %d, Nz: %d, Npca: %d, nR: %d, nP: %d, HJ: %d\n',...
    N,Ny,Nx,Nz,Npca,nR,nP,HuhJhun);
alpha = 0.05;
fprintf('Results:\n')

fprintf('FPR (parametric):\n');      disp(mean(ppararep<=alpha));
fprintf('FPR (permutation [F]):\n'); disp(mean(ppermrep<=alpha));
fprintf('FPR (permutation [r]):\n'); disp(mean(pcorrrep<=alpha));
fprintf('Mean CCs (unpermuted)\n');  disp(mean(corrFirst));
fprintf('Mean CCs (random perm)\n'); disp(mean(corrLast));

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

% =================================================================
function Q = BLUSmtx(X)
% Create (N-P) x N matrix that generates BLUS residuals

N = size(X,1);
P = size(X,2);

% Find safe subset to drop
Done=false;
while (~Done)
  I1 = randperm(N,N-P);
  I0 = setdiff(1:N,I1);
  if rank(X(I0,:))==P
    Done=true;
  end
end

% Magnus, J. R., & Sinha, A. K. (2005). On Theil’s errors. The Econometrics
% Journal, 8(1), 39–54.

% S and Q (A) are transposed relative to Magnus & Sinha
M    = eye(N) - X*inv(X'*X)*X';
S    = eye(N);
S    = S(I1,:);
Q    = sqrtm(inv(S*M*S'))*S*M;

