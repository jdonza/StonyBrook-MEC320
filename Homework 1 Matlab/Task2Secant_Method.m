close all
clear
clc
%Secant Method
%Equation
x = .1;
h = .1;
fx = (20*exp(-x/5))*cos(x)*sin(x)*(x-7) + erfc(x) - 20;
xi = x + h;
fxi = (20*exp(-xi/5))*cos(xi)*sin(xi)*(xi-7) + erfc(xi) - 20;

%Incremental Search
while fx*fxi > 0
    fx = (20*exp(-x/5))*cos(x)*sin(x)*(x-7) + erfc(x) - 20;
    h = .1;
    xi = x + h;
    fxi = (20*exp(-xi/5))*cos(xi)*sin(xi)*(xi-7) + erfc(xi) - 20;
    x = xi;
end

%When x = 1.8, fx = -3.9348
%When x = 1.9, fx = 1.3469
%Initial Guesses
xjj = 1.8; % xj-1
xj= 1.9; % xj
Es = .001; %Stopping Criterion
Ea = abs(((xj - xjj)/(xj))*100);
i = 0;
while Ea > Es
    i = i+1;
    fj = (20*exp(-xj/5))*cos(xj)*sin(xj)*(xj-7) + erfc(xj) - 20;
    fjj = (20*exp(-xjj/5))*cos(xjj)*sin(xjj)*(xjj-7) + erfc(xjj) - 20;
    xjjj = xj - ((fj*(xj-xjj))/(fj - fjj)); %Definition of secant (xj+1)
    Ea = abs(((xjjj - xj)/(xjjj))*100);
    xj = xjjj;
end
disp('Root 1 =');
disp(xjjj);


%Start at initial guess of 2 and do the same thing

%Equation
x = 2;
h = .1;
fx = (20*exp(-x/5))*cos(x)*sin(x)*(x-7) + erfc(x) - 20;
xi = x + h;
fxi = (20*exp(-xi/5))*cos(xi)*sin(xi)*(xi-7) + erfc(xi) - 20;

%Incremental Search
while fx*fxi > 0
    fx = (20*exp(-x/5))*cos(x)*sin(x)*(x-7) + erfc(x) - 20;
    h = .1;
    xi = x + h;
    fxi = (20*exp(-xi/5))*cos(xi)*sin(xi)*(xi-7) + erfc(xi) - 20;
    x = xi;
end

%When x = 2.6, fx = 3.1104
%When x = 2.7, fx = -0.6358
%Initial Guesses
xjj = 2.6; % xj-1
xj= 2.7; % xj
Es = .001; %Stopping Criterion
Ea = abs(((xj - xjj)/(xj))*100);
ii = 0;
while Ea > Es
    ii = ii+1;
    fj = (20*exp(-xj/5))*cos(xj)*sin(xj)*(xj-7) + erfc(xj) - 20;
    fjj = (20*exp(-xjj/5))*cos(xjj)*sin(xjj)*(xjj-7) + erfc(xjj) - 20;
    xjjj = xj - ((fj*(xj-xjj))/(fj - fjj)); %Definition of secant (xj+1)
    Ea = abs(((xjjj - xj)/(xjjj))*100);
    xj = xjjj;
end
disp('Root 2 =');
disp(xjjj);


