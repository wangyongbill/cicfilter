clear
clc
close all
%% 非2的整数次幂
T = 1;
fs = 240;
fm1 = 100;
fm2 = 110;
fft_N = 400;
N = fft_N*2.56;

t = 0:1/fs:(N-1)/fs;
Len = length(t);
y = sin(2*pi*fm1*t)+sin(2*pi*fm2*t);

Y = abs(fft(y));
Y(1) = Y(1)/Len;
Y(2:Len) = 2*Y(2:Len)/Len;

f = (0:Len-1)*fs/Len;
figure;plot(f,Y)

%% 补零到2的整数次幂
n = 2^nextpow2(Len);
Y = abs(fft(y,n));
Y(1) = Y(1)/n;
Y(2:n) = 2*Y(2:n)/n;

f = (0:n-1)*fs/n;
figure;plot(f(1:n/2+1),Y(1:n/2+1))