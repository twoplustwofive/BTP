function [U,D, fronorm, spec_norm] = singlepassRSVD_Hermitian(A,k,p)
    [m,n] = size(A);
    G = randn(n, k+p);
    Y = A*G;
    [Q,~] = qr(Y);
    C = (Q*Y)/(Q*G);
    [Uhat, D] = eig(C);
    U = Q*Uhat;
    fronorm = norm(A-Q*Q'*A,"fro");
    spec_norm = norm(A-Q*Q'*A);
end

