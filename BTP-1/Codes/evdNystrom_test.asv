clear all;
clc;

k = 100;
p = 10;
A = generateSPDmatrix(1000);

[U1, lambda] = evdNystrom(A,k,p);

[U2, D2, V2] = svd(A);
% rhs = (1+(k/(p-1)));
% 
% singular_square_root = 0;
% 
% for i = k+1:1000
%     singular_square_root = singular_square_root+D2(i,i)*D2(i,i);
% end
% rhs = rhs*(singular_square_root);
% rhs = sqrt(rhs);
% disp(fronorm+" "+rhs);

