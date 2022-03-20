close all
clear
clc

x0 = 0;
y0 = 0;
%Parabolic Interpolation
%Equations
f = @(x,y)10*sin(x+3) - 6*cos((2*y) + 2);
df_dx = @(x,y)10*cos(x+3);
df_dy = @(x,y)12*sin((2*y) + 2);
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
    b = g(h1);
    c = g(h2);
    
    h3 = (a * (h1^2 - h2^2) +  b * (h2^2 - x0^2) + c * (h0^2 - h1^2)) / (2*a*(h1-h2) + 2*b*(h2-h0) + 2*c*(h0-h1));
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








