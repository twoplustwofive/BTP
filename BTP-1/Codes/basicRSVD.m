function [U, D, V, fronorm, spec_norm] = basicRSVD(A,k,p)
    [m,n] = size(A);
    G = randn(n, k+p);
    Y = A*G;
    [Q,~] = qr(Y,0);
    B = Q'*A;
    [Uhat, D, V] = svd(B, 'econ');
    U = Q*Uhat;
    fronorm = norm(A-Q*Q'*A,"fro");
    spec_norm = norm(A-Q*Q'*A);
end

