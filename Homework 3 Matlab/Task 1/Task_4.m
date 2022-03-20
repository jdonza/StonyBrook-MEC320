close all
clear
clc

%%%%%%%%%% Trapezoidal Rule

%Step 1: Extract the pressure and volume data from excel

P_kpa = xlsread('Data.xlsx',4,'E3:E3602'); %Extracts pressure data in kpa
V_m3 = xlsread('Data.xlsx',4,'D3:D3602'); %Extracts volume data in m^3

%Step 2: Convert pressure data into pascals

[n,m] = size(P_kpa);
P_pa = zeros(n,m);
for i = 1:n
    P_pa(i,m) = P_kpa(i,m)/1000;
end

%Step 3: Calculate the work in Joules over the whole cycle

Work_j = zeros(n,m);
for i = 1:n-1
    Work_j(i,m) = (V_m3(i+1,m) - V_m3(i,m)) * ((P_pa(i,m) + P_pa(i+1,m))/2); %Calculates the work in joules
end
S = sum(Work_j);
%Step 4: Plot the pressure as a function of volume in matlab

figure (1)
plot(V_m3,P_pa); %Plot made with linear x and y axes
title('Pressure vs. Volume (Linear)');
xlabel('Volume m^3')
ylabel('Pressure Pa')
grid on
figure (2)
loglog(V_m3,P_pa); %Plot made with logarithmic axes
title('Pressure vs. Volume (Logarithmic)');
xlabel('Volume m^3')
ylabel('Pressure Pa')

grid on

%Step 5: Plot the volume versus the crank angle

crank_angle = xlsread('Data.xlsx',4,'C3:C3602'); %Grabs crank angle data from excel
figure (3)
plot(crank_angle,V_m3);
title ('Volume vs. Crank Angle');
xlabel('Crank Angle in Degrees');
ylabel('Volume m^3');
grid on

%Step 6: Use a central finite difference to find the derivative of the
%volume

V_prime = zeros(n,m);
h = .2; %Defines step size of crank angle data
for i = 2:n-1
    V_prime(i,m) = (V_m3(i+1,m) - V_m3(i-1,m)) / (2*h);
end
figure (4)
plot (crank_angle,V_prime);
title('Derivative of Volume vs. Crank Angle');
xlabel('Crank Angle')
ylabel('Derivative of Volume')
grid on

%Step 7: Implement an incremental search to estimate the range in which the
%zeros are located

vi = 1;
vii = 1;
j = 0;
for x = 2:n-2
    vi= V_prime(x,m);
    vii = V_prime(x+2,m);
    if vi*vii < 0
        j = j+1;
        lower(j,m) = crank_angle(x);
        lower(j,m+1) = x; %Puts cell number of value in next column
        upper(j,m) = crank_angle(x+2);
        upper(j,m+1) = x+2; %Puts cell number of value in next column
    end
end

%Step 8: Determine an approximate sine function to model the graph

max_v = max(V_prime);
min_v = min(V_prime);
max_crank = max(crank_angle);
min_crank = min(crank_angle);

for x = 1:n
    y(x,m) = max_v*sind(crank_angle(x,m));
end
figure (5)
plot(crank_angle,y);
title('Approximate Sine Function to Represent Derivative of Volume')
xlabel('Crank Angle')
ylabel('Derivative of Volume')
grid on

fprintf('The sine function to approximately fit the data is y = %f*sinx\n',max_v)

%Step 9: Utilizing the approximate sine function and the bisection method,
%determine the location of the root.

ea = 100;
es = .0005; %units in percent
i = 0;
for k = 1:j
    while ea > es
        i = i+1;
        xu = upper(k,m);
        xl = lower(k,m);
        Xr = (xu+xl)/2;
        func_Xr = max_v*(sind(Xr));
        func_xl = max_v*(sind(xl));
        func_xu = max_v*(sind(xu));
        
        if func_xl*func_Xr < 0 %If true, the root lies in the lower sub-interval
            xu = Xr;
        elseif func_xl*func_Xr > 0 %If true, the root lies in the upper sub-interval
            xl = Xr;
        else %Only true if func_xl*func_Xr = 0
            xroot(k,m) = Xr;
            break
        end
        ea = abs((xu-xl)/(xu+xl)) *100;
    end
end
disp('Crank Angle in which root occurs =');
disp(xroot);

%Step 10: Calculate the work for each individual stroke

for i = 1:lower(i,m+1)
    k = 1;
    Work = zeros(lower(i,m+1),m);
    while V_prime(k,m) >= 0
        Work(k,m) = (V_m3(k+1,m) - V_m3(k,m)) * ((P_pa(k,m) + P_pa(k+1,m))/2); %Calculates the work in joules
        k = k+1;
    end
    while V_prime(k,m) <=0
        Work(k,m) = (V_m3(k+1,m) - V_m3(k,m)) * ((P_pa(k,m) + P_pa(k+1,m))/2); %Calculates the work in joules
        k = k+1;
    end 
end





