function [Z,T,bh,sig2] = glm_miss(X,Y,M,c)
    % Compute GLM allowing for missingness in X or Y, for multiple Y
    % X - NxP design
    % Y - NxK data
    % M - NxK mask - 1 if present, 0 if missing
    % c - 1xP contrast
    %
    % Depends on mtimesx package 
    % https://uk.mathworks.com/matlabcentral/fileexchange/25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support
    % No NaN's are allowed, rather, the presence missingness in a given column of Y 
    % must be indicated in the mask matrix M

    global MX MY XtMX XtMXi % Keep these to improve memory management with permutation

    if any(isnan(Y(:)))
        error('NaNs not allowed; encode missingness with M')
    end

    N = size(X,1);
    P = size(X,2);
    K = size(Y,2);
    
    nM = sum(M,1);

    % Mask model, data
    MY = reshape(M.*Y,[N,1,K]);
    MX = reshape(M,[N,1,K]).*X;

    % Compute OLS estimates accounting for missingness
    XtMX = mtimesx(X',MX);
    XtMXi = zeros(size(XtMX));
    for k = 1:K
        XtMXi(:,:,k) = inv(XtMX(:,:,k));
    end
    XtMY = mtimesx(X',MY);
    bh = mtimesx(XtMXi,XtMY);
    DF  = nM - P;

    sig2 = sum((MY - mtimesx(MX,bh)).^2,1) ./ reshape(DF,[1 1 K]);
    
    con = mtimesx(c,bh);
    SE2 = mtimesx(c,mtimesx(XtMXi,c')).*sig2;
    T   = con./sqrt(SE2);
    T   = reshape(T,[1,K]);

    Z   = zeros(1,K);
    % Upper tail t to z
    I = Z>=0;
    Z(I) = -norminv(tcdf(T(I),DF(I),'upper'));
    % Lower tail t to z
    I = Z<0;
    Z(I) = norminv(tcdf(T(I),DF(I)));

return
