function [U, D, V, fronorm, spec_norm] = singlepassRSVD(A,k,p)
    [m,n] = size(A);
    Gc = randn(n,k+p);
    Gr = randn(m,k+p);
    Yc = A*Gc;
    Yr = A*Gr;
    [Qc,~] = qr(Yc);
    [Qr,~] = qr(Yr);
%     disp(size(Gr)+" "+size(Qc)+" "+size(Yr)+" "+size(Qr)+" "+size(Qc)+" "+size(Yc)+" "+size(Qr)+" "+size(Gc));
    C = lsqminnorm(Gr'*Qc, Yr'*Qr);
%     disp(size(C));
    [Uhat, D, Vhat] = svd(C);
    U = Qc*Uhat;
    V = Qr*Vhat;
    Q = Qc;
    fronorm = norm(A-Q*Q'*A,"fro");
    spec_norm = norm(A-Q*Q'*A);
end

