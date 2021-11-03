clear all;
clc;

x_k = [];
y_fronorm = [];
y_spectral = [];
y_theo_spectral = [];
y_theo_fronorm = [];
im_num = 1;
for k = 100:5:500
    p = 10;
    A = randi([-1,1],1000,1000);
    [U1, D1, V1, fronorm, spec_norm] = basicRSVD(A,k,p);
    [U2, D2, V2] = svd(A);
    rhs = (1+(k/(p-1)));
    
    singular_square_root = 0;
    
    for i = k+1:1000
        singular_square_root = singular_square_root+D2(i,i)*D2(i,i);
    end
    rhs = rhs*(singular_square_root);
    rhs = sqrt(rhs);
%     disp(fronorm+" "+rhs);
    
    % average error in spectral norm
    rhs_spectral = 0;
    rhs_spectral = (1 + sqrt(k/(p-1)));
    for i = k+1:1000
        if i == k + 1
            reqd_sigma = D2(i,i)*D2(i*i);
        end
        singular_square_root = singular_square_root+D2(i,i)*D2(i,i);
    end
    
    rhs_spectral = rhs_spectral*reqd_sigma;
    rhs_spectral = rhs_spectral + ((2.17*sqrt(k+p))/p)*sqrt(singular_square_root);
%     rhs_spectral = sprintf('%0.7f',rhs_spectral);
%     disp(spec_norm + " " + rhs_spectral);
    x_k(end+1) = k;
    y_spectral(end+1) = spec_norm;
    y_fronorm(end+1) = fronorm;
    y_theo_fronorm(end+1) = rhs;
    y_theo_spectral(end+1) = rhs_spectral;
end
figure;
plot(x_k, y_fronorm);
xlabel('Target rank - k');
ylabel('Frobenius norm');
title('Target rank Vs Error');
% saveas(gcf, sprintf('plots/basicRSVD_test_%d.png', im_num)); im_num = im_num + 1;


figure;
plot(x_k, y_spectral);
xlabel('Target rank - k');
ylabel('Spectral norm');
title('Target rank Vs Error');
% saveas(gcf, sprintf('plots/basicRSVD_test_%d.png', im_num)); im_num = im_num + 1;


disp("Average Observed frobrenius norm = "+mean(y_fronorm));
disp("Average Theoretical frobrenius norm = "+mean(y_theo_fronorm));
disp("Average Observed Spectral norm = "+mean(y_spectral));
disp("Average Theoretical Spectral norm = "+mean(y_theo_spectral));

x_k = [];
y_fronorm = [];
y_spectral = [];
y_theo_spectral = [];
y_theo_fronorm = [];
k = 500;
for n = 1000:10:1500
    p = 10;
    A = randi([-1,1],n,n);
    [U1, D1, V1, fronorm, spec_norm] = basicRSVD(A,k,p);
    [U2, D2, V2] = svd(A);
    rhs = (1+(k/(p-1)));
    
    singular_square_root = 0;
    
    for i = k+1:n
        singular_square_root = singular_square_root+D2(i,i)*D2(i,i);
    end
    rhs = rhs*(singular_square_root);
    rhs = sqrt(rhs);
%     disp(fronorm+" "+rhs);
    
    % average error in spectral norm
    rhs_spectral = (1 + sqrt(k/(p-1)));
    for i = k+1:n
        if i == k + 1
            reqd_sigma = D2(i,i)*D2(i*i);
        end
        singular_square_root = singular_square_root+D2(i,i)*D2(i,i);
    end
    
    rhs_spectral = rhs_spectral*reqd_sigma;
    rhs_spectral = rhs_spectral + ((2.17*sqrt(k+p))/p)*sqrt(singular_square_root);
%     rhs_spectral = sprintf('%0.7f',rhs_spectral);
%     disp(spec_norm + " " + rhs_spectral);
    x_k(end+1) = n;
    y_spectral(end+1) = spec_norm;
    y_fronorm(end+1) = fronorm;
    y_theo_fronorm(end+1) = rhs;
    y_theo_spectral(end+1) = rhs_spectral;
end

figure;
plot(x_k, y_fronorm);
xlabel('Matrix size(n)');
ylabel('Frobenius norm');
title('Matrix size(n) Vs Error');
% saveas(gcf, sprintf('plots/basicRSVD_test_%d.png', im_num)); im_num = im_num + 1;

figure;
plot(x_k, y_spectral);
xlabel('Matrix size(n)');
ylabel('Spectral norm');
title('Matrix size(n) Vs Error');
% saveas(gcf, sprintf('plots/basicRSVD_test_%d.png', im_num)); im_num = im_num + 1;

disp("Average Observed frobrenius norm = "+mean(y_fronorm));
disp("Average Theoretical frobrenius norm = "+mean(y_theo_fronorm));
disp("Average Observed Spectral norm = "+mean(y_spectral));
disp("Average Theoretical Spectral norm = "+mean(y_theo_fronorm));


