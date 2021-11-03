function [U,lambda] = evdNystrom(A,k,p)
    [m,n] = size(A);
    G = randn(n, k+p);
    Y = A*G;
    [Q,~] = qr(Y);
    B1 = A*Q;
    B2 = Q'*B1;
    C = chol(B2);
    F = B1/C;
    [U,sigma,~] = svd(F,'econ');
    lambda = sigma*sigma;
end

