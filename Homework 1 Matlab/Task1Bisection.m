close all
clear
clc
%Bisection Method
%Initial Lower and Upper Bound Guess

xu = 7;
xl = 1;

%Begin loop for 8 iterations
i = 0;
for n = 1:8
    i = i+1;
    Xr = (xu+xl)/2;
    func_Xr = (sqrt(Xr)*sin(Xr))-1;
    func_xl = (sqrt(xl)*sin(xl))-1;
    func_xu = (sqrt(xu)*sin(xu))-1;
    
    if func_xl*func_Xr < 0 %If true, the root lies in the lower sub-interval
        xu = Xr;
    elseif func_xl*func_Xr > 0 %If true, the root lies in the upper sub-interval
        xl = Xr;
    else %Only true if func_xl*func_Xr = 0
        xroot = Xr;
    end
    ea = abs((xu-xl)/(xu+xl)) *100;
end
disp('Final Root =');
disp(Xr);
disp('Approximate Error =')
disp(ea)
