%%%TASK 1
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

close all
clear
clc

%%Gauss Elimination with partial pivoting

%The general concepts of how to construct partial pivoting in Matlab was obtained from 
%the pseudocode provided in the textbook on page 266. This part of the code was interpreted 
%into Matlab and altered to better represent this particular task.

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

%Step 5: Foward Elimination with partial pivoting

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

%Step 6: Back Substitution

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
        

%%%TASK 2

close all
clear
clc

%%LU-Decomposition with iterative refinement

%Information pertaining to the implementation of the decomposition of the A
%matrix was obtained from the pseudocode presented on page 282 in the
%textbook.

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

%Step 2: Determine the size of the matrix

[nA,mA] = size(A);%nA = rows; mA = columns
[nb,mb] = size(b);%nb = rows; mb = columns

%Step 3: Begin Decomposing A into L and U matrix

L = zeros(nA,mA);
U = zeros(nA,mA);

%Step 3.1: Foward Elimination to Determine U Matrix
for j = 1:nb %j is for columns
    for i = 1:nA
        U(i,j) = A(i,j);
    end
end

for j = 1:nb
    for i = j+1:nA %i is for row
        U(i,:) =U(i,:)-(U(j,:)* U(i,j)/U(j,j));
    end
end

%Step 3.2: Determine L Matrix

for i = 1:nA
    L(i,i) = 1;
end

for j = 1:nb
    for i = j+1:nA
        L(i,j) = A(i,j)/A(j,j);
    end
    
end

%Step 4: Begin substitution to find solution matrix

%Step 4.1: Foward substitution with L matrix

d = zeros(n,1); %Sets up a zero matrix to fill in solution values
m = 1; 
d(m,:) = b(m,:)/L(m,m);
for i = m+1:1:nA
    solution = b(i,:);
    for j = 1:nb
        solution = solution - L(i,j)*d(j);
    end
    d(i) = solution/L(i,i);
end
disp('[d_approx] =');
disp(d);

%Step 4.2: Backwards substitution with U matrix

x_approx = zeros(n,1); %Sets up a zero matrix to fill in approximate solution values
x_approx(n) = d(n,:)/U(n,n); %Defines the value for your last variable in matrix
for i = n-1:-1:1
    solution = d(i,:); %Calls upon solution value as given in b-vector
    for j = i+1:n
        solution = solution - U(i,j)*x_approx(j,:);
    end
    x_approx(i) = solution/U(i,i);
end
disp('[x] =');
disp(x_approx);

%Step 5: Iterative Refinement

es = zeros(n,1);
es(:,1) = .001; %units in percent
eA = zeros(n,1); 
eA(:,1) = 100; %units in percent
k = 0;

while max(eA) > es(:,1)
    k = k+1;
    b_approx = A * x_approx;
    disp('b_approx =');
    disp(b_approx);
    E = b - b_approx; %Error in b matrix
    
    %Find the correction factor matrix x_delta
    
    %Foward substitution with L matrix
    
    y = zeros(n,1); %Sets up a zero matrix to fill in solution values
    m = 1;
    y(m,:) = E(m,:)/L(m,m);
    for i = m+1:1:nA
        solution = E(i,:);
        for j = 1:nb
            solution = solution - L(i,j)*y(j);
        end
        y(i) = solution/L(i,i);
    end
    disp('[y] =');
    disp(y);
    
    %Backwards substitution with U matrix
    
    x_delta = zeros(n,1); %Sets up a zero matrix to fill in approximate solution values
    x_delta(n) = y(n,:)/U(n,n); %Defines the value for your last variable in matrix
    for i = n-1:-1:1
        solution = y(i,:); %Calls upon solution value as given in b-vector
        for j = i+1:n
            solution = solution - U(i,j)*x_delta(j,:);
        end
        x_delta(i) = solution/U(i,i);
    end
    disp('[x_delta] =');
    disp(x_delta);
    
    x_approx1 = x_delta + x_approx;
    eA = abs((x_approx1 - x_approx)./(x_approx1))*100;
    x_approx = x_approx1;
end
x = x_approx;
disp('x =');
disp(x);



%%%TASK 3

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



%%%TASK 4

tic
close all
clear
clc

%%Steepest Ascent/Descent
%Define anonymous functions
x0 = 0;
y0 = 0;
Ea = zeros(2,1);
Ea(:,:) = 100; %units in percent
I = max(Ea);
Es = .01; %units in percent
w = 0;

    while I > Es
        
        w = w+1;
        f = @(x,y)10*sin(x+3) - 6*cos(2*y + 2);
        df_dx = @(x,y)10*cos(x+3);
        df_dy = @(x,y)12*sin(2*y + 2);
        df2_dx2 = @(x,y)-10*sin(x + 3);
        df2_dy2 = @(x,y)24*cos(2*y + 2);
        df2_dxdy2 = @(x,y)10*cos(x+3) + 12*sin(2*y + 2);
        
        g = @(h) f(x0 + (df_dx(x0,y0))*h, y0 + (df_dy(x0,y0))*h);
        
        %Initial Guesses
        h0 = 0;
        h1 = 0.1;
        h2 = 0.2;
        
        ea = 100;
        es = .01; %units in percent
        k = 0;
        while ea > es
            k = k+1;
            a = g(h0);
            I = g(h1);
            c = g(h2);
            
            h3 = (a * (h1^2 - h2^2) +  I * (h2^2 - x0^2) + c * (h0^2 - h1^2)) / (2*a*(h1-h2) + 2*I*(h2-h0) + 2*c*(h0-h1));
            h_vector(k,:) = h3;
            d = g(h3);
            h0 = h1;
            h1 = h2;
            h2 = h3;
            
            if k > 1
                ea = abs((h_vector(k,:) - h_vector(k-1,:))/(h_vector(k,:)))*100;
            end
        end
        disp('h_opt =');
        disp (h3);
        
        x0_old = x0;
        y0_old = y0;
        x0 = x0 + df_dx(x0,y0)*h3;
        y0 = y0 + df_dy(x0,y0)*h3;
        
        Ea1 = abs((x0 - x0_old)/(x0))*100;
        Ea2 = abs((y0 - y0_old)/(y0))*100;
        Ea(1,1) = Ea1;
        Ea(2,1) = Ea2;
        I = max(Ea);
        
    end
    disp(x0);
    disp(y0);
    toc
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    




