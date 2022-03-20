close all
clear
clc

%Incremental Search
i = 0;
for x = 0:1:20 
    fx = (50 * exp(-x/5) * cos(x)) - 2;
    h = 1;
    xi = x + h;
    fxi = (50 * exp(-xi/5) * cos(xi)) - 2;
    if fx*fxi < 0
        i = i+1;
        upper(i) = xi; %This stores my upper bound limits into an array
        lower(i) = x; %This stores my lower bound limits into an array
    end
end


%Begin Root Finding Methods
%Define Stopping Criterion
Es = 5; %Units in Percent
for r = 1:i
    
    lowerbound = lower(r);
    upperbound = upper(r);
    Ea = 100; %No error for first iteration, put value at 100 so loop will perform
    j=0;
    
    %Begin False Position Method 
    
    while Ea > Es
        j = j+1;
        fxl1 = (50 * exp(-lowerbound/5) * cos(lowerbound)) - 2;
        fxu1 = (50 * exp(-upperbound/5) * cos(upperbound)) - 2;
        xr = upperbound - ((fxu1*(lowerbound-upperbound))/(fxl1 - fxu1));
        fxr = (50 * exp(-xr/5) * cos(xr)) - 2;
        
        if fxl1 * fxr > 0
            lowerbound = xr; %Root lies in upper sub-interval
        elseif fxl1 * fxr < 0
            upperbound = xr; %Root lies in lower sub-interval
        else %Only true if fxl1*fxr = 0
            xroot = xr;
        end
        Ea = abs((upperbound - lowerbound)/(upperbound + lowerbound))*100;      
    end
    disp('Number of Iterations for False Position Method =');
    disp(j);
    
    %Being Modified Secant Method
    
    Esnew = .00005; %Units in percent
    x0 = xr;
    delta = 0.01;
    Eanew = 100;
    k = 0;
    while Eanew > Esnew
        k = k+1;
        fx0 = (50 * exp(-x0/5) * cos(x0)) - 2;
        x1 = x0 + x0*delta;
        fx1 = (50 * exp(-x1/5) * cos(x1)) - 2;
        x11 = x0 - ((delta*fx0*x0)/(fx1-fx0));
        Eanew = abs((x11 - x0)/(x11))*100;
        x0 = x11;
    end 
    disp('Number of Iterations for the Modified Secant Method =');
    disp(k);
    disp('Root =')
    disp(x11);
    
end 
    