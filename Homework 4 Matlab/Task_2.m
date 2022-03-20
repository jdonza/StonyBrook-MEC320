close all
clear
clc

%%%%%%%%%%Euler's Method and Bisection Method (Part A)

%Step 1: Define Boundary Conditions

To = 10; %Units of degrees celsius
Tf = 275; %T(17) in units of degress celsius
h_prime2 = 5.5*10^-8;
Ta = 30; %Units of degrees celsius
h = 0.001; %Step size for Euler's Method
Es = 0.01; %Units in percent (Stopping Criterion)
Ea = 100; %Units in percent (Initial approximate error)
iter = 17/h; %Number of iterations needed to evaluate over entire interval
x0 = 0; %Starting point on x-axis
zl = 6; %Lower initial guess for Bisection
zu = 8; %Upper initial guess for Bisection
k = 0; %Counter for Bisection Iterations

%Step 2: Define anonymous functions for Euler's Method

Ti = @(T,z) T + z*h; %Equation 1
zi = @(T,z) z + (h_prime2*(T-Ta)^4)*(h); %Equation 2

%Step 3: Create Loop for Euler's Method

while Ea > Es
    %Perform Bisection Method
    zr = (zl + zu)/2;
    %Perform Euler's Method
    for i = 1:iter
        if i == 1
            x_a(i,1) = x0;
        else
            x_a(i,1) = x_a(i-1,1) + h;
        end
        if i == 1
            T_a(i) = (Ti(To,zr));
            z_a(i) = (zi(To,zr));
        else
            T_a(i) = (Ti(T_a(i-1),z_a(i-1)));
            z_a(i) = (zi(T_a(i-1),z_a(i-1)));
        end
    end
    %Continue with Bisection Method
    Ea = abs((Tf - T_a(end))/(Tf))*100;
    g = Tf - T_a(end);
    if g > 0
        zl = zr;
    elseif g < 0
        zu = zr;
    else
        zroot = zr;
    end
    k = k+1;
end
figure (1)
plot (x_a,T_a);
title 'Temperature vs. X-Values'
xlabel 'X-Values'
ylabel 'Temperature (Celsius)'
hold on
%%%%%%%%%%Begin Part B
zl = 6;
zu = 8;
Ea = 100;
Es = .01;
hnew = input('Coarser Step Size =')
iter1 = 17/hnew;
p = 0;
%Define new anonymous functions
Tii = @(T,z) T + z*hnew; %Equation 1
zii = @(T,z) z + (h_prime2*(T-Ta)^4)*(hnew); %Equation 2
while Ea > Es
    %Perform Bisection Method
    zr = (zl + zu)/2;
    %Perform Euler's Method
    for i = 1:iter1
        if i == 1
            x_b(i,1) = x0;
        else
            x_b(i,1) = x_b(i-1,1) + hnew;
        end
        if i == 1
            T_b(i) = (Tii(To,zr));
            z_b(i) = (zii(To,zr));
        else
            T_b(i) = (Tii(T_b(i-1),z_b(i-1)));
            z_b(i) = (zii(T_b(i-1),z_b(i-1)));
        end
    end
    %Continue with Bisection Method
    Ea = abs((Tf - T_b(end))/(Tf))*100;
    g = Tf - T_b(end);
    if g > 0
        zl = zr;
    elseif g < 0
        zu = zr;
    else
        zroot = zr;
    end
    p = p+1;
end
[n,m] = size(T_b);
[x_val,indx] = intersect(roundn(x_a,-5),roundn(x_b,-5));
T_d = T_a(indx);
for s = 1:m
    E(s) = abs((T_d(s) - T_b(s))/(T_d(s)))*100;
end
w = max(E);
disp('Maximum Error =');
disp(w);
if w > 0.5
    disp('Use smaller step size')
elseif w < 0.5
    disp('Use greater step size')
else
    disp('Step size achieved error of 0.5% exactly')
end

figure (1)
plot (x_b,T_b);
%%%%%%%%%%Begin Part C

%Step 2: Define New Constants



%Step 1: Define K-functions
Ea_c = 100; %Units in percent
Es_c = .01; %Units in percent
z(1) = 6;

while Ea_c > Es_c
    for i = 1:iter
        if i == 1
            x_c(i,1) = x0;
        else
            x_c(i,1) = x_c(i-1,1) + h;
        end
        if i == 1
            k1_T(i) = Ti(To,z);%K11
            k1_z(i) = zi(To,z);%K12
            
            a = To + ((k1_T(i))*(h/2));
            b = z + ((k1_z(i))*(h/2));
            
            k2_T(i) = Ti(a,b);
            k2_z(i) = zi(a,b);
            
            c = To + ((k2_T(i))*(h/2));
            d = z + ((k2_z(i))*(h/2));
            
            k3_T(i) = Ti(c,d);
            k3_z(i) = zi(c,d);
            
            e = To + ((k3_T(i))*h);
            f = z + ((k3_z(i))*h);
            
            k4_T(i) = Ti(e,f);
            k4_z(i) = zi(e,f);
            
            T_c (i) = To + ((1/6) * (k1_T(i) + 2*k2_T(i) + 2*k3_T(i) + k4_T(i)))*h;
            z_c(i) = z + ((1/6) * (k1_z(i) + (2*k2_z(i)) + (2*k3_z(i)) +k4_z(i)))*h;
        else
            k1_T(i) = Ti(T_c(i-1),z_c(i-1));%K11
            k1_z(i) = zi(T_c(i-1),z_c(i-1));%K12
            
            a(i) = T_c(i-1) + ((k1_T(i))*(h/2));
            b(i) = z_c(i-1) + ((k1_z(i))*(h/2));
            
            k2_T(i) = Ti(a(i),b(i));
            k2_z(i) = zi(a(i),b(i));
            
            c(i) = T_c(i-1) + ((k2_T(i))*(h/2));
            d(i) = z_c(i-1) + ((k2_z(i))*(h/2));
            
            k3_T(i) = Ti(c(i),d(i));
            k3_z(i) = zi(c(i),d(i));
            
            e(i) = T_c(i-1) + ((k3_T(i))*h);
            f(i) = z_c(i-1) + ((k3_z(i))*h);
            
            k4_T(i) = Ti(e(i),f(i));
            k4_z(i) = zi(e(i),f(i));
            
            T_c (i) = T_c(i-1) + ((1/6) * (k1_T(i) + (2*k2_T(i)) + (2*k3_T(i)) +k4_T(i)))*h;
            z_c(i) = z_c(i-1) + ((1/6) * (k1_z(i) + (2*k2_z(i)) + (2*k3_z(i)) +k4_z(i)))*h;
            
        end
    end
    Ea_c(i) = abs((Tf - T_c(end))/(Tf))*100;
    g = Tf - T_c(end);
end

Ea_c = 100; %Units in percent
Es_c = .01; %Units in percent
z(1) = 8;

while Ea_c > Es_c
    for i = 1:iter
        if i == 1
            x_c(i,1) = x0;
        else
            x_c(i,1) = x_c(i-1,1) + h;
        end
        if i == 1
            k1_T(i) = Ti(To,z);%K11
            k1_z(i) = zi(To,z);%K12
            
            a = To + ((k1_T(i))*(h/2));
            b = z + ((k1_z(i))*(h/2));
            
            k2_T(i) = Ti(a,b);
            k2_z(i) = zi(a,b);
            
            c = To + ((k2_T(i))*(h/2));
            d = z + ((k2_z(i))*(h/2));
            
            k3_T(i) = Ti(c,d);
            k3_z(i) = zi(c,d);
            
            e = To + ((k3_T(i))*h);
            f = z + ((k3_z(i))*h);
            
            k4_T(i) = Ti(e,f);
            k4_z(i) = zi(e,f);
            
            T_c (i) = To + ((1/6) * (k1_T(i) + 2*k2_T(i) + 2*k3_T(i) + k4_T(i)))*h;
            z_c(i) = z + ((1/6) * (k1_z(i) + (2*k2_z(i)) + (2*k3_z(i)) +k4_z(i)))*h;
        else
            k1_T(i) = Ti(T_c(i-1),z_c(i-1));%K11
            k1_z(i) = zi(T_c(i-1),z_c(i-1));%K12
            
            a(i) = T_c(i-1) + ((k1_T(i))*(h/2));
            b(i) = z_c(i-1) + ((k1_z(i))*(h/2));
            
            k2_T(i) = Ti(a(i),b(i));
            k2_z(i) = zi(a(i),b(i));
            
            c(i) = T_c(i-1) + ((k2_T(i))*(h/2));
            d(i) = z_c(i-1) + ((k2_z(i))*(h/2));
            
            k3_T(i) = Ti(c(i),d(i));
            k3_z(i) = zi(c(i),d(i));
            
            e(i) = T_c(i-1) + ((k3_T(i))*h);
            f(i) = z_c(i-1) + ((k3_z(i))*h);
            
            k4_T(i) = Ti(e(i),f(i));
            k4_z(i) = zi(e(i),f(i));
            
            T_c (i) = T_c(i-1) + ((1/6) * (k1_T(i) + (2*k2_T(i)) + (2*k3_T(i)) +k4_T(i)))*h;
            z_c(i) = z_c(i-1) + ((1/6) * (k1_z(i) + (2*k2_z(i)) + (2*k3_z(i)) +k4_z(i)))*h;
            
        end
    end
    Ea_c(i) = abs((Tf - T_c(end))/(Tf))*100;
    g = Tf - T_c(end);
end


