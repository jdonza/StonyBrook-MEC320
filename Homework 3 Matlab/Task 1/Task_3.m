close all
clc
clear

%%%%%%%%%% Inverse Fourier Transform

%Step 1: Extract data from sheet one in excel file into matlab matrix

[num,A] = xlsread('Data.xlsx',3);
Fk = str2double(A);
[n,m] = size(Fk);

%Step 2: Calculate the power spectrum

for i = 1:n
    pk(i,:) = 2*abs(Fk(i))^2;
end
% The sampling frequency of the audio file is 44,100 Hz. The data is
% equi-spaced in the frequency domain starting at zero and ending at the
% sampling frequency. Therefore, the frequency is incremented at
% 44,100/34,032 = 1.2958 Hz per pk.
increment = 44100/n;

for i = 1:n
    if i == 1
        xi(i,:) = 0;
    else
        xi(i,:) = xi(i-1,:) + increment;
        if xi(i,:) > 800
            break
        end
    end
end
figure(1)
plot (xi,pk(1:i))
xlabel('0 Hz to 800 Hz');
ylabel('Decibel')
title('Power Spectrum')
grid on
i = i+1;
for j = i:n
    xi(j,:) = xi(j-1,:) + increment;
    if xi(j,:) > 3500
        break
    end
end
figure(2)
plot (xi(i:j),pk(i:j))
xlabel('0 Hz to 3500 Hz');
ylabel('Decibels');
title('Power Spectrum')
grid on

% Take the inverse fourier transform
T = 1/44100;
w0 = (2*pi)/T;

for i = 1:n-1
    if i == 1
        xj(i,:) = 0;
    else
        xj(i,:) = xj(i-1,:) + increment;
    end
end



for g = 1:n-1
    if g == 1
        ti(g,:) = 0;
    elseif g == 2
        ti(g,:) = 1/increment;
    else
        ti(g,:) = 1/(xj(g-1,:)) + (1/increment);
    end 
end 

for k = 1:n-1
    fn(k,1) = 0;
    %imaginary(k,1) = 0;
    for s = 0:n-1
        angle = k*w0*s;
        fn(k,1) = fn(k,1) + (Fk(k,1)*cos(angle)) - (1i*Fk(k,1)*sin(angle)) ;
        %imaginary(k,1) = imaginary(k,1) - (1i*Fk(k,1)*sin(angle));
    end
end

figure (3)
plot(ti,fn)
title('Amplitude Data vs. Time')
grid on


  
        
