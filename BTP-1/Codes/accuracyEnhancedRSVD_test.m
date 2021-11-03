clear all;
clc;

k = 100;
p = 10;
A = randi([-1,1],1000,1000);
[m,n] = size(A);

[U1, D1, V1, specform] = accuracyEnhancedRSVD(A,k,p,1);

[U2, D2, V2] = svd(A);

rhs = (1 + sqrt(k/(p-1))) + ((2.71828*sqrt(k+p))/p)*(sqrt(min(n,m) - k));
rhs = rhs.^(1/(2*1 + 1));
% reqd_sigma = 0;
i = k + 1;
reqd_sigma = D2(i,i)*D2(i,i);
% end
rhs = rhs * (reqd_sigma);
% rhs = sqrt(rhs);
% fronorm = sprintf('%0.7f',fronorm);
disp(specform+" "+rhs);




x_k = [];
y_spectral = [];
y_theo_spectral = [];
im_num = 1;
for k = 100:5:500
    p = 10;
    A = randi([-1,1],1000,1000);
    [U1, D1, V1, specnorm] = accuracyEnhancedRSVD(A,k,p,1);
    [U2, D2, V2] = svd(A);
    rhs = (1 + sqrt(k/(p-1))) + ((2.71828*sqrt(k+p))/p)*(sqrt(min(n,m) - k));
    rhs = rhs.^(1/(2*1 + 1));
    % reqd_sigma = 0;
    i = k + 1;
    reqd_sigma = D2(i,i)*D2(i,i);
    % end
    rhs = rhs * (reqd_sigma);
    % rhs = sqrt(rhs);
    % fronorm = sprintf('%0.7f',fronorm);
%     disp(specform+" "+rhs);
    x_k(end+1) = k;
    y_spectral(end+1) = specnorm;
    y_theo_spectral(end+1) = rhs;
end

figure;
plot(x_k, y_spectral);
xlabel('Target rank - k');
ylabel('Spectral norm');
title('Target rank Vs Error');
% saveas(gcf, sprintf('plots/basicRSVD_test_%d.png', im_num)); im_num = im_num + 1;

disp("Average Observed Spectral norm = "+mean(y_spectral));
disp("Average Theoretical Spectral norm = "+mean(y_theo_spectral));

x_k = [];
y_spectral = [];
y_theo_spectral = [];
k = 500;
for n = 1000:10:1500
    p = 10;
    A = randi([-1,1],n,n);
    [U1, D1, V1, specnorm] = accuracyEnhancedRSVD(A,k,p,1);
    [U2, D2, V2] = svd(A);
    rhs = (1 + sqrt(k/(p-1))) + ((2.71828*sqrt(k+p))/p)*(sqrt(min(n,m) - k));
    rhs = rhs.^(1/(2*1 + 1));
    % reqd_sigma = 0;
    i = k + 1;
    reqd_sigma = D2(i,i)*D2(i,i);
    % end
    rhs = rhs * (reqd_sigma);
    % rhs = sqrt(rhs);
    % fronorm = sprintf('%0.7f',fronorm);
%     disp(specform+" "+rhs);
    x_k(end+1) = n;
    y_spectral(end+1) = specnorm;
    y_theo_spectral(end+1) = rhs;
end


figure;
plot(x_k, y_spectral);
xlabel('Matrix size(n)');
ylabel('Spectral norm');
title('Matrix size(n) Vs Error');
% saveas(gcf, sprintf('plots/basicRSVD_test_%d.png', im_num)); im_num = im_num + 1;

disp("Average Observed Spectral norm = "+mean(y_spectral));
disp("Average Theoretical Spectral norm = "+mean(y_theo_spectral));
