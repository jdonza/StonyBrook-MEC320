x = linspace(0,10);
y = (20.*exp(-x/5)).*cos(x).*sin(x).*(x-7) + erfc(x) - 20;
plot(x,y)
grid on