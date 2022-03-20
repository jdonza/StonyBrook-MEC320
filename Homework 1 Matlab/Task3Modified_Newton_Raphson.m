close all
clear
clc
%Modified Newton Raphson Method
%Initial Guess
x0 = 0.1;

%Define stopping criterion
Es = 0.001;  %Units in percent

%Define Equations
xi = x0;
fx = xi^2 + 6 - exp(2*xi);
fxprime = 2*xi - 2*exp(2*xi); % First derivative of function
fxdprime = 2 - 4*exp(2*xi); % Second derivative of function
xii = xi - (fx*fxprime)/((fxprime^2) - (fx*fxdprime));
Ea = abs((xii - xi)/(xii))*100;
i = 0;
while Ea > Es
    fx = xi^2 + 6 - exp(2*xi);
    fxprime = 2*xi - 2*exp(2*xi); % First derivative of function
    fxdprime = 2 - 4*exp(2*xi); % Second derivative of function
    xii = xi - (fx*fxprime)/((fxprime^2) - (fx*fxdprime));
    Ea = abs((xii - xi)/(xii))*100;
    xi = xii;
    i = i+1;
end
disp('Root =');
disp(xii);


