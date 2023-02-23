clear
clc
close all

Fs = 160e6; % sample rate
R = 4000;  % decimator factor
D = 1;  % differential delay
N = 1;  % number of stage
Fp = 7.5e3; % pass band
Fstp = 15e3; % stop band
Ap = 0.01; % attenuation in pass band
Astp = 60; % attenuation in stop band

CICDecim = dsp.CICDecimator(R, D, N);
CICCompDecim = dsp.CICCompensationDecimator(CICDecim, ...
    'DecimationFactor',R,'PassbandFrequency',Fp, ...
    'StopbandFrequency',Fstp,'SampleRate',Fs/R);
fvtool(CICDecim,CICCompDecim,...
cascade(CICDecim,CICCompDecim),'ShowReference','off','Fs',[Fs Fs/R Fs])
legend('CIC Decimator','CIC Compensator','Resulting Cascade Filter');


