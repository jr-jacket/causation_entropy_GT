function [px] = p_mkde(x,X,h)
%This code performs KDE. Takes as input the point of interest (x), data
%used for estimation (X), and the bandwidth to be used (h). This file
%should not be changed unless replacing the density estimation method. 


[d N] = size(X);
covX = cov(X');
detS = abs(det(covX));
invS = pinv(covX);
sum = 0;
for i = 1:N
    K = (x - X(:,i))'*invS*(x - X(:,i));
    sum = sum + 1/(sqrt((2*pi)^d*detS)*h^d)*exp(-K/(2*h^2));
end
px = 1/(N)*sum;
return