function [U,D,V, specnorm] = accuracyEnhancedRSVD(A,k,p,q)
    [m,n] = size(A);
    G = randn(n, k+p);
    Y = A*G;
    for j = 1:q
        Z = A'*Y;
        Y = A*Z;
    end
    [Q, R] = qr(Y,0);  %modified; earlier command Q = qr(A) was incorrect
    B = Q'*A;
    [Uhat, D, V] = svd(B, 'econ');
    U = Q*Uhat;
    specnorm = norm(A-Q*Q'*A);
end

