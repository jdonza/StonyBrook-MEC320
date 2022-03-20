close all
clear 
clc

%Step 1: Define Anonymous Function

fx = @(x) (x^0.1) * (1.18-x) * (1-(exp(20*(x-1))));

%Step 2: Define Uppper and Lower x bounds

x0 = 0.05;
xf = 1.2;

%Step 3: Define intial h1 and h2 values

h1 = 0.05;
h2 = 0.025;
h = [h1;h2]; %Creates step size vector

%Step 4: Define Error Criterion

Es = 0.1; %Units in percent

%Step 5: Calculate O(h^2) using centered finite difference
k = 1;
j = 0;
Ea = 100;
f = 0; %Sets counter for plot for part B
x= x0;
p = 0;
while x < xf
h = [h1;h2]; %Creates step size vector
p = p+1;
i = 0; %iteration number
Ea = 100;
D(1,1) = (fx(x+h(1,1)) - fx(x-h(1,1)))/(2*h(1,1));
while Ea > 0.1 
    i= i+1;
    D(i+1,1) = (fx(x+h(i+1,1)) - fx(x-h(i+1,1)))/(2*h(i+1,1)); %CFD
 
        for k = 2:i+1
            j = 2 + i - k;
            
            D(j,k) = (((4^(k-1)) * (D(j+1,k-1))) - (D(j,k-1))) / (4^(k-1) - 1); %Richardson's Extrapolation
        end
        Ea(j,1) = abs((D(j,k) - D(j+1,k-1))/ (D(j,k)))*100; %Calculates Percent Error
        if Ea > Es
            g = size(h,1);
            h(g+1,1) = h(g,1)/2;
        end
end
f = f+1;
iter(p) = i; 
derv(p) = x; %Counts the x-values (Points where derivative is calculated at) 
xii(f) = h(k,1) + derv(p); %Calculates x-value one step size greater than x-value used to calculate derivative
xi(f) = derv(p) - h(k,1); %Calculates x-value one step size smaller than x-value used to calculate derivative
di(f) = D(1,k);
g = size(h,1);
x = x+h(g,1);
end

%Plots for part A
figure (1)
plot(derv,iter);
title 'Number of Iterations vs. X-Values'
xlabel 'X-Values'
ylabel 'Number of Iterations'
grid on

%Plots for part B
[n,m] = size(derv);
for i = 1:m
    func(1,i) = fx(derv(1,i));
    func1(1,i) = fx(xi(1,i)); %Evaluates function at xi-1
    func2(1,i) = fx(xii(1,i)); %Evalutes function at xi+1
end 
figure (2)
plot(derv,func,'-*','MarkerFaceColor','blue');
title 'Plot of Function'
xlabel 'X-values'
ylabel 'Function Value'
grid on

hold on
scatter(derv,func1,'.','red');

hold on
scatter(derv,func2,'.','red');
legend('Points Used to Calculate','Points Used to Evaluate');

%Plot for part D
figure (3)
plot(derv,di);
title 'Derivative vs. X-Values'
xlabel 'X-Values'
ylabel 'Derivative of Function'
grid on
