close all
clear
clc

%%Gauss Seidel Method

%Step 1: Create Matrix from system of equations

A = [2 -1 -2 9 2 1;...
    2 2 -2 3 11 1;...
    1 1 9 -2 3 1;...
    13 3 -4 -2 1 2;...
    -2 -1 1 4 2 11;...
    -4 -12 1 2 2 1];

b = [7; -12; 16; 13; 4; 8];

n = input('Number of variables =');

%Step 2: Define Size of the matrix

[nA,mA] = size(A);%nA = rows; mA = columns
[nb,mb] = size(b);%nb = rows; mb = columns

%Step 3: Check for convergence

for i = 1:mA
        k = abs(A(i,i));
        s = sum(abs(A(i,:))- abs(k)); %Checking if A-Matrix is diagnolly dominant
        if s > k
            disp('System may diverge');
            break
        end 
end 


%Shows that intially, the system will not converge
% Insert Diagonally Dominant Matrix

A = [13 3 -4 -2 1 2; ...
    -4 -12 1 2 2 1; ...
    1 1 9 -2 3 1;...
    2 -1 -2 9 2 1; ...
    2 2 -2 3 11 1; ...
    -2 -1 1 4 2 11];

%Step 4: Define your initial guesses

x_initial = zeros(mA,1);
x = zeros(mA,1);

%Step 5: Begin Gauss-Seidel Iteration

es = zeros(nb,1);
es(:,:) = 0.0005; %Units in percent
ea = zeros(nb,1);
ea(:,:) = 100;
k = 0;

while  max(ea) > es(:,1)
    k = k+1; %Calculates number of iterations
    for i = 1:nA
        solution = b(i,:);
        for j = 1:i-1
            solution = solution - A(i,j)*x(j,:);
        end
        for j = i+1:mA
            solution = solution - A(i,j)*x(j,:);
        end
        x(i) = solution/A(i,i);
    end
    ea(:,1) = abs((x(:,1) - x_initial(:,1))./(x(:,1)))*100;
    x_initial = x;
end
disp('x =');
disp(x);



