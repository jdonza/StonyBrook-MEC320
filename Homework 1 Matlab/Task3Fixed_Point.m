close all
clear
clc
%Fixed Point Method
%Initial Guess
x0 = 0.1;

%Define stopping criterion
Es = 0.001;  %Units in percent
xi = (log(x0^2+6))/2;
xii = (log(xi^2+6))/2;
gprime = abs((xi - xii)/(x0 - xi));
Ea = abs((xii - xi)/(xii))*100;
i = 0;
while Ea > Es
    i = i+1;
    xi = (log(x0^2+6))/2;
    xii = (log(xi^2+6))/2;
    Ea(i) = abs((xii - xi)/(xii))*100;
    gprime = abs((xi - xii)/(x0 - xi))
    x0 = xii;
end
disp('g(x) is linearly converging because the absolute value of gprime is less than one.');
disp('Root =');
disp(xii);
semilogy(Ea);

