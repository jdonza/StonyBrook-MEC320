close all
clear
clc

%%LU-Decomposition with iterative refinement

%Information pertaining to the implementation of the decomposition of the A
%matrix was obtained from the pseudocode presented on page 282 in the
%textbook.

%Step 1: Create Matrix from system of equations
A = [2 -5 1; -1 3 -1; 3 -4 2];

b = [12; -8; 16];

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






