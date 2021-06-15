function [Z,T,bh,sig2] = glm_miss(X,Y,M,c)
    % Compute GLM allowing for missingness in X or Y, for multiple Y
    % X - NxP design
    % Y - NxK data
    % M - NxK mask - 1 if present, 0 if missing
    % c - 1xP contrast
    %
    % No NaN's are allowed, rather, the presence any missingness in Y 
    % must be indicated in the mask matrix M.  (Missingness in X could 
    % be accommodated, by indicating missingness in any row i of X with
    % *all* K elements of row i of M being set to zero, but has no
    % advantage over just removing that row from the outset.)
    % 
    % If pagemtimes isn't available, depends on mtimesx package 
    % https://uk.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support
    %________________________________________________________________________
    % TE Nichols June 2021

    Use_pinv=0;

    if exist('pagemtimes')
        mymtimes    = @(x,y)pagemtimes(x,y);
        mymtimest1  = @(x,y)pagemtimes(x,'transpose',y,'none');
        mymtimest2  = @(x,y)pagemtimes(x,'none',y,'transpose');
    elseif exist('mtimesx')
        mymtimes    = @(x,y)mtimesx(x,y);
        mymtimest1  = @(x,y)mtimesx(x,'t',y);
        mymtimest2  = @(x,y)mtimesx(x,y,'t');
    else
        error('Neither pagemtimes or mtimesx exists; upgrade to Matlab R2020b or install mtimesx')
    end

    % Make big variables permanent to improve memory management with permutation
    global MX MY pMX MXtMXi bh sig2 con SE2 T

    if any(isnan(Y(:)))
        error('NaNs not allowed; encode missingness with M')
    end

    N = size(X,1);
    P = size(X,2);
    K = size(Y,2);
    
    % Mask model, data
    MY = reshape(M.*Y,[N,1,K]);              % MY:     N x 1 x K
    MX = reshape(M,[N,1,K]).*X;              % MX:     N x P x K

    % Check if contrast is estimable for every k
    for k = 1:K
        I  = M(:,k)~=0;
        mX = squeeze(MX(I,:,k));
        if any(max(abs(c'-mX'*pinv(mX')*c'))>1e-8)
            warning(sprintf('Missingness in element %d leads to unestimable contrast',k))
        end
    end

    % Compute OLS estimates accounting for missingness, 
    if Use_pinv
        % "betahat = pinv(MX)*MY"
        pMX = zeros(P,N,K);                 % pMX:     P x N x K
        for k = 1:K
            I = M(:,k)~=0;
            pMX(:,I,k) = pinv(squeeze(MX(I,:,k)));
        end
        bh = mymtimes(pMX,MY);              % bh:      P x 1 x K
    else
        % "betahat = inv(MX'MX)*MX'*MY"
        MXtMXi = zeros(P,P,K)               % MXtMXi:  P x N x K
        for k=1:K
            MXtMXi(:,:,k) = inv(MX(:,:,k)'*MX(:,:,k));
        end
        bh = mymtimes(MXtMXi,mymtimest1(MX,'t',MY));
                                            % bh:      P x 1 x K
    end

    % Inference...               sig2,con,SE2,T,Z:   1 x 1 x K
    DF   = sum(M)-P;
    sig2 = sum((MY - mymtimes(MX,bh)).^2,1) ./ reshape(DF,[1 1 K]);
    con  = mymtimes(c,bh);
    SE2  = mymtimes(c,mymtimes(mymtimest2(pMX,pMX),c')).*sig2;
    T    = con./sqrt(SE2);
    T    = reshape(T,[1,K]);

    Z = zeros(1,K);
    % Upper tail t to z
    I = Z>=0;
    Z(I) = -norminv(tcdf(T(I),DF(I),'upper'));
    % Lower tail t to z
    I = Z<0;
    Z(I) = norminv(tcdf(T(I),DF(I)));

return
