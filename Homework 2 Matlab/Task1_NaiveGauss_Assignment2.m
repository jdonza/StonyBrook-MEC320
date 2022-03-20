close all
clear
clc

%%Gauss Elimination without partial pivoting

%The general concepts of how to set-up the for loops for the foward elimination
%and back substitution was obtained from the pseudocode provided in the
%textbook on page 254. This part of the code was interpreted into Matlab and altered to
%better represent this particular task.

%Step 1: Create Matrix from system of equations
A = [-2 1 -9 2 2 5 1 8 -1; ...
    1 -6 3 12 11 -3 -7 2 -2; ...
    2 -3 -1 2 5 13 2 16 3; ...
    -6 2 -3 1 2 -1 14 2 -2; ...
    -1 -13 4 -4 3 7 -3 2 -3; ...
    15 3 -5 7 -2 5 2 1 -1; ...
    14 3 -5 7 -2 5 2 1 -1; ...
    3 2 2 -5 -2 -14 6 7 -2; ...
    3 2 3 -4 -2 -13 6 7 -2];

b = [6; 10; 25; -12; 18; 23; 24; 8; 9];

n = input('Number of variables =');

%Step 2: Create augmented matrix

Am = [A b];

%Step 3: Determine the size of the matrix

[nA,mA] = size(A);%nA = rows; mA = columns
[nb,mb] = size(b);%nb = rows; mb = columns

%Step 4: Determine the condition number

S = abs(sum(A)); %Finds the absolute sum of all the columns
S_max = max(S); %Determines the maximum sum 
I = inv(A); %Calculates the inverse of A
I_sum = abs(sum(I)); %Sums the columns in I
I_max = max(I_sum); %Chooses the maximum sum
cond_A = S_max * I_max; %Determines the condition number

if cond_A > 100
    disp('system is ill-conditioned');
else 
    disp('system is well-conditioned');
end 

%Step 5: Foward Elimination 

for j = 1:nb %j is for columns
    for i = j+1:nA %i is for rows
        Am(i,:) =Am(i,:)-(Am(j,:)* Am(i,j)/Am(j,j));
    end
end

% Step 6: Back Substitution

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

%Step 7: Check to see if code works

b = A*x;
disp('b =');
disp(b);
