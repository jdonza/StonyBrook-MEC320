close all
clear
clc

%Initial Lower and Upper Bound Guess

xu = 7;
xl = 1;

for n = 1:8
    Xr = (xu+xl)/2;
    func_Xr = (sqrt(Xr)*sin(Xr))-1;
    func_xl = (sqrt(xl)*sin(xl))-1;
    func_xu = (sqrt(xu)*sin(xu))-1;
    
    if func_xl*func_Xr < 0
        xu = Xr;
    elseif func_xl*func_Xr > 0
        xl = Xr;
    else
        xroot = Xr;
    end
    
end
