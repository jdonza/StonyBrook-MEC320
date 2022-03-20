close all
clear
clc
%Equation
%Incremental Search
i = 0;
for x = 0:.1:20 %Interval from 0-2
    fx = (50 * exp(-x/5) * cos(x)) - 2;
    h = .1;
    xi = x + h;
    fxi = (50 * exp(-xi/5) * cos(xi)) - 2;
    if fx*fxi < 0
            i = i+1;
        upper(i) = xi;
        lower(i) = x;
    end
end
