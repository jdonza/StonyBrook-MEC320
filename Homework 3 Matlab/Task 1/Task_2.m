close all
clear
clc

%%%%%%%%%%Linear Spline
%Step 1: Extract data from sheet one in excel file into matlab matrix

A = xlsread('Data.xlsx',2);

%Step 2: Determine the size of matrix A

[n,u] = size(A);

%Step 3: Place values from A matric in matrix f and x

for i = 1:n
    f(i,:) = A(i,2);
    x(i,:) = A(i,1);
end


%Step 4: Determine the function values at each x-value utilizing the slopes determined in
%step 3.

for i = 1:n-1
    r = linspace(x(i),x(i+1),1000);
    m(i,:) = (f(i+1) - f(i)) / (x(i+1) - x(i));
    fx = @(x) f(i,:) + m(i,:)*(x-x(i));
    fprintf('The linear spline equation to fit the data is y = %f + %f*x\n',f(i,:),m(i,:))
    plot(x,f,'-o',r,fx(r))
    title('Linear Spline')
    xlabel('x-values')
    ylabel('y-values')
    hold on
end

