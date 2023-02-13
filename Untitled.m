close all
clear all
clc
 
%set system parameter
fs = 2500;    %The frequency of the local oscillator signal
Fs = 160e6;   %sampling frequency
N =  24;         %Quantitative bits
L = 5000;
 
%Generating an input signal
t =0:1/Fs:(1/Fs)*(L-1);          %Generating the time series of sampling frequencies
sc = sin(2*pi*fs*t);        %a sinusoidal input signal that produces a random starting phase
 
b =[1,-1];%integerator
a =[1,-1];%comb 
 
%comb
c1=filter(1,b,sc);
c2=filter(1,b,c1);
c3=filter(1,b,c2);
 
y = downsample(c3,4);
 
%integerater
i1 =filter(a,1,y);
i2 =filter(a,1,i1);
i3 =filter(a,1,i2);
sf = i3;
 
% N = 3 R = 4 M = 1
coe = [1 1 1 1];
coe1 = conv(coe,coe);
coe2 = conv(coe1,coe);
sf1 =filter(coe2,1,sc);
sf2 = downsample(sf1,4);
 
figure(1),
subplot(211),stem(t(1:1600),sc(1:1600));
xlabel('时间(t)','fontsize',8);
ylabel('幅度(dB)','fontsize',8);
title('sc','fontsize',8);
subplot(212),stem(t(1:400),sf(1:400));
xlabel('时间/4(t)','fontsize',8);
ylabel('幅度(dB)','fontsize',8);
title('sf','fontsize',8);

figure(2),
subplot(211),stem(t(1:1600),sc(1:1600));
xlabel('时间(t)','fontsize',8);
ylabel('幅度(dB)','fontsize',8);
title('sc','fontsize',8);
subplot(212),stem(t(1:400),sf2(1:400));
xlabel('时间/4(t)','fontsize',8);
ylabel('幅度(dB)','fontsize',8);
title('sf','fontsize',8);

figure(3)
plot(sf2 - sf) 
figure(3)
plot(sf2 - sf)

