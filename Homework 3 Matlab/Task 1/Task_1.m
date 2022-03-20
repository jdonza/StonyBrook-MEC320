close all
clear
clc

%%%%%%%%%%Linear Least-Squares Regression
%Step 1: Extract data from sheet one in excel file into matlab matrix

A = xlsread('Data.xlsx',1);

%Step 2: Define unknowns in a1 and a0 equations

[u,m] = size(A);
sum_x = 0;
sum_y = 0;
sum_xy = 0;
sum_x2 = 0;
for i = 1:u
    sum_x = sum_x + A(i,1);
    sum_y = sum_y + A(i,2);
    sum_xy = sum_xy + (A(i,1)*A(i,2));
    sum_x2 = sum_x2 + (A(i,1)*A(i,1));
end

%Step 3: Solve for a1 and a0

a1 = ((u*sum_xy) - (sum_x*sum_y)) / ((u*sum_x2) - (sum_x)^2);
y_average = sum_y/u;
x_average = sum_x/u;
a0 = y_average - (a1*x_average);

%Step 4: Display the linear regressed equation

fprintf('The linear regressed equation to fit the data is y = %f + %f*x\n',a0,a1)

%Step 5: Find correlation coefficient

St = 0;
Sr = 0;
for i = 1:u
    St = St + (A(i,2)-y_average)^2;
    Sr = Sr + (A(i,2)-a0-(a1*A(i,1)))^2;
end

r2 = (St-Sr)/(St); %Coefficient of determination
r = sqrt(r2); %Correlation coefficient
disp('Correlation Coefficient for Linear Regressed Model = ');
disp(r);
disp('Coefficient of determination for the Linear Regressed Model =');
disp(r2);

%Step 6: Determine the new y-values for the linear regressed model 

for i = 1:u
    C(i,1) = a0 + a1*A(i,1);
end 

%%%%%%%%%%Parabolic Least-Squares Regression

%Step 1: Define unknowns for normal equations

sum_x = 0;
sum_y = 0;
sum_xy = 0;
sum_x2 = 0;
sum_x3 = 0;
sum_x4 = 0;
sum_x2y = 0;

for i = 1:u
    sum_x = sum_x + A(i,1);
    sum_y = sum_y + A(i,2);
    sum_xy = sum_xy + (A(i,1)*A(i,2));
    sum_x2 = sum_x2 + (A(i,1)*A(i,1));
    sum_x3 = sum_x3 + (A(i,1))^3;
    sum_x4 = sum_x4 + (A(i,1))^4;
    sum_x2y = sum_x2y + ((A(i,1))^2)*A(i,2);
end

%Step 2: Use symbolic toolbox to set up matrix. x = a0, y = a1, z = a2.

syms x y z
eqn1 = (u*x) + (sum_x*y) + (sum_x2*z) == sum_y;
eqn2 = (sum_x*x) + (sum_x2*y) + (sum_x3*z) == sum_xy;
eqn3 = (sum_x2*x) + (sum_x3*y) + (sum_x4*z) == sum_x2y;
[a,b] = equationsToMatrix([eqn1, eqn2, eqn3],[x, y, z]);

%Step 3: Use Gauss Elimination with partial pivoting to solve for a0, a1, and a2

n = input('Number of unknowns for Gauss Elimination =');%In this particular case, n=3.

%Step 4: Create augmented matrix

Am = [a b];

%Step 5: Determine the size of the matrix

[nA,mA] = size(a);%nA = rows; mA = columns
[nb,mb] = size(b);%nb = rows; mb = columns


%Step 6: Foward Elimination with partial pivoting

for j = 1:nb %j is for columns
    %Partial Pivoting
    m = n+1;
    p = 1;
    k = p;
    pivot = Am(k,k);
    for ii = k+1:1:n
        pivot2 = Am(ii,k); 
        if pivot2 > pivot
            pivot = pivot2;
            p = ii;
        end
    end
    if p ~= k
        for jj = k:1:n
            pivot2 = Am(p,jj);
            Am(p,jj) = Am(k,jj);
            Am(k,jj) = pivot2;
        end
        pivot2 = Am(p,m);
        Am(p,m) = Am(k,m);
        Am(k,m) = pivot2;
    end
    
    for i = j+1:nA %i is for rows
        Am(i,:) =Am(i,:)-(Am(j,:)* Am(i,j)/Am(j,j));
    end
end

%Step 7: Back Substitution

x = zeros(n,1); %Sets up a zero matrix to fill in solution values
m = n+1; %Written to insure we call last column in augmented matrix 
x(n) = Am(n,m)/Am(n,n); %Defines the value for your last variable in matrix
for i = n-1:-1:1 
   solution = Am(i,m); %Calls upon solution value as given in b-vector
   for j = i+1:n 
       solution = solution - Am(i,j)*x(j,:);
   end 
   x(i) = solution/Am(i,i);
end
disp('[x] =');
disp(x);

%Step 8: Check to see if code works

b = a*x;
disp('b =');
disp(vpa(b));        

%Step 9: Display Parabolic Equation

fprintf('The polynomial regressed equation to fit the data is y = %f + %f*x + %f*x^2\n',x(1),x(2),x(3))

%Step 10: Find correlation coefficient

St1 = 0;
Sr1 = 0;
for i = 1:n
    St1 = St1 + (A(i,2)-y_average)^2;
    Sr1 = Sr1 + (A(i,2) - x(1) - (x(2)*A(i,1)) - (x(3)*(A(i,2))^2))^2;
end

r2_1 = (St1-Sr1)/(St1); %Coefficient of determination
r_1 = sqrt(r2_1); %Correlation coefficient
disp('Correlation Coefficient for the Polynomial Regressed Model = ');
disp(r_1);
disp('Coefficient of determination for the Polynomial Regressed Model =');
disp(r2_1);

%Step 11: Determine the new y-values for the parabola

for i = 1:u
    D(i,1) = x(1,1) + x(2,1)*A(i,1) + (x(3,1)*(A(i,1))^2);
end 

%%%%%%%%%% Custom Sine Regression Model

%Step 1: Linearize the sine function and data into the form of
%sin^-1(y) = A + Bx where A = a0 and B = a1

for i = 1:u
    ylin(i,:) = asin(A(i,2));
    xlin(i,:) = A(i,1);
end 
F = [xlin ylin]; %Puts the data into one matrix for simplicity

%Step 2: Define unknowns in a1 and a0 equations

[u,m] = size(F);
sum_x = 0;
sum_y = 0;
sum_xy = 0;
sum_x2 = 0;
for i = 1:u
    sum_x = sum_x + F(i,1);
    sum_y = sum_y + F(i,2);
    sum_xy = sum_xy + (F(i,1)*F(i,2));
    sum_x2 = sum_x2 + (F(i,1)*F(i,1));
end

%Step 3: Solve for a1 and a0

a1 = ((u*sum_xy) - (sum_x*sum_y)) / ((u*sum_x2) - (sum_x)^2);
y_average = sum_y/u;
x_average = sum_x/u;
a0 = y_average - (a1*x_average);

%Step 4: Display the custom sine equation

fprintf('The custom sine function to fit the data is y = sin(%f + %f*x)\n',a0,a1)

%Step 5: Find correlation coefficient

St = 0;
Sr = 0;
for i = 1:u
    St = St + (A(i,2)-y_average)^2;
    Sr = Sr + (A(i,2)-a0-(a1*A(i,1)))^2;
end

r2 = (St-Sr)/(St); %Coefficient of determination
r = sqrt(r2); %Correlation coefficient
disp('Correlation Coefficient for Custom Sine fit = ');
disp(r);
disp('Coefficient of determination for the Custom Sine fit =');
disp(r2);

%Step 6: Determine the new y-values for the linear regressed model 

for i = 1:u
    G(i,1) = sin(a0 + (a1 * A(i,1)));
end 

%%%%%%%%%%Exponential Least Sqaures Regression

%Step 1: Linearize the exponential function and data into the form of
%log(y) = log(alpha) + beta*x; where log(alpha)=a0, beta=a1,
%log(y)=ylin

for i = 1:u
    ylin(i,:) = log(A(i,2));
    xlin(i,:) = A(i,1);
end 
B = [xlin ylin]; %Puts the data into one matrix for simplicity

%Step 2: Define unknowns in a1 and a0 equations

[u,m] = size(B);
sum_x = 0;
sum_y = 0;
sum_xy = 0;
sum_x2 = 0;
for i = 1:u
    sum_x = sum_x + B(i,1);
    sum_y = sum_y + B(i,2);
    sum_xy = sum_xy + (B(i,1)*B(i,2));
    sum_x2 = sum_x2 + (B(i,1)*B(i,1));
end

%Step 3: Solve for a1 and a0

a1 = ((u*sum_xy) - (sum_x*sum_y)) / ((u*sum_x2) - (sum_x)^2);
y_average = sum_y/u;
x_average = sum_x/u;
a0 = y_average - (a1*x_average);

%Step 4: Solve for alpha and beta

beta = a1;
alpha = exp(a0);

%Step 5: Display the Exponential regressed equation

fprintf('The exponential regressed equation to fit the data is y = %f*e^%f*x\n',alpha,beta)

%Step 6: Find correlation coefficient

St = 0;
Sr = 0;
for i = 1:u
    St = St + (A(i,2)-y_average)^2;
    Sr = Sr + (A(i,2)-a0-(a1*A(i,1)))^2;
end

r2 = (St-Sr)/(St); %Coefficient of determination
r = sqrt(r2); %Correlation coefficient
disp('Correlation Coefficient for the Exponential Model = ');
disp(r);
disp('Coefficient of determination for the Exponential Model =')
disp(r2);

%Step 7:Determine the new y-values for the linear regressed model 

for i = 1:u
    E(i,1) = alpha*exp(beta*A(i,1));
end 


%%%%%%%%%% Plot all the data
figure
subplot(2,2,1)
plot(A(:,1),A(:,2),'x',A(:,1),C(:,1));
grid on
xlabel('x-values')
ylabel('y-values')
title('Linear Regression Plot')

subplot(2,2,2)
plot(A(:,1),A(:,2),'x',A(:,1),D(:,1));
grid on
xlabel('x-values')
ylabel('y-values')
title('Parabolic Regression Plot');

subplot(2,2,3)
plot(A(:,1),A(:,2),'x',A(:,1),E(:,1));
grid on
xlabel('x-values')
ylabel('y-values')
title('Exponential Regression Plot');

subplot(2,2,4)
plot(A(:,1),A(:,2),'x',A(:,1),G(:,1));
grid on
xlabel('x-values')
ylabel('y-values')
title('Custom Sine Regression Plot');


