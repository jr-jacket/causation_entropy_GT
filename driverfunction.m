function driverfunction 
close all
[xdata,time]=makedyn; %Generate discrete dynamics

potentialFunctions=@(x,t)[x(1,:); x(2,:);sin(5*t);x(1,:).^2]; %anonymous function defining the function vector F considered by the CEM
%Every row is a potential function with each column representing a
%subsequent data point. 

perm=100; %number of permutations to be used in the permutation test

[rawCEM,adjustedCEM,Tmatrix]=calccem(xdata,time,potentialFunctions,perm);   %Compute CEM
rawCEM
adjustedCEM             %display raw and adjusted CEM 

return

function [x,timevect]=makedyn %Function to make dynamics according to system in readme from randomized initial conditions for 350 data points
numbpoints=350;
A=-[1.5 -1.25;.95 -.35];
T=0.01;
IC=randn(2,1);
x=zeros(2,numbpoints);
x(:,1)=IC;

timevect=zeros(1,length(x));
for i=2:numbpoints
    timevect(i)=timevect(i-1)+T;
     x(:,i)=T*(A*x(:,i-1)+[0;sin(5*timevect(i-1))])+x(:,i-1);      
end
T*A+eye(2)
size(x)
plot(timevect,x)
xlabel('Time (s)')
ylabel('State Value')
legend('$x_1$','$x_2$', 'interpreter','latex')


return