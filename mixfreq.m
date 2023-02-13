%% 混频后的滤波器
clear
clc
close all
%% 采集设置
fs = 160e6;
ts = 1/fs;
T = 0.015;
t = 0:ts:T-ts;
len = length(t);

fpass = 8e6; 
fstop = 10e6;
Num = GetLowFilterCoef(fs,fpass,fstop,0.01,100,20);

%% 未抽取的噪声
noise = randn(1,length(t));
% noise = diff(noise);
power = sum(noise*noise')/length(noise);
figure;plot(noise);title( ['原始噪声波形,功率为:' num2str(power)])

[noise_abs,f_noise] = Func_pufenxi(noise,fs,0);
figure;plot(f_noise,noise_abs);title(['原始噪声频谱,功率为:' num2str(power)])

%% 直接抽取的噪声
lag = 8;
noise_chouqu = noise(1:lag:end);
figure;plot(noise_chouqu);
power = sum(noise_chouqu*noise_chouqu')/length(noise_chouqu);
title( ['未滤波8抽取的噪声波形,功率为:' num2str(power)])

[noise_abs_chouqu,f_noise_chouqu] = Func_pufenxi(noise_chouqu,fs/lag,0);
figure;plot(f_noise_chouqu,noise_abs_chouqu);
title( ['未滤波8抽取的噪声频谱,功率为:' num2str(power)])

%% 求和平均后抽取
fpass = 15e6; 
fstop = fpass+5e6;
Num = GetLowFilterCoef(fs,fpass,fstop,0.1,40,20);
% figure;fvtool(Num,1)
% noise = filter(Num,1,noise);

sumNum = 80;
noise_Mat = reshape(noise,[sumNum,length(noise)/sumNum]);
noise_sumNum = sum(noise_Mat,1)/sumNum;

power = sum(noise_sumNum*noise_sumNum')/length(noise_sumNum);

figure;
plot(noise_sumNum);title( ['8求和平均抽取的噪声波形,功率为:' num2str(power)])

[noise_abs_sum,f_noise_sum] = Func_pufenxi(noise_sumNum,fs/sumNum,0);
figure;plot(f_noise_sum,noise_abs_sum);
title( ['8求和平均抽取的噪声频谱,功率为:' num2str(power)])

%% 噪声经过低通滤波器后
noise = filter(Num,1,noise);
power = sum(noise*noise')/length(noise);

figure;plot(noise);
title( ['滤波后噪声波形,功率为:' num2str(power)])

[noise_abs,f_noise] = Func_pufenxi(noise,fs,0);
figure;plot(f_noise,noise_abs);
title( ['滤波后噪声频谱,功率为:' num2str(power)])

%% 噪声经低通滤波后抽取
lag = 8;
noise_chouqu = noise(1:lag:end);
power = sum(noise_chouqu*noise_chouqu')/length(noise_chouqu);
figure;plot(noise_chouqu);
title( ['滤波后8抽取的噪声波形,功率为:' num2str(power)])

[noise_abs_chouqu,f_noise_chouqu] = Func_pufenxi(noise_chouqu,fs/lag,0);
figure;plot(f_noise_chouqu,noise_abs_chouqu);
title( ['滤波后8抽取的噪声频谱,功率为:' num2str(power)])

%% CIC滤波器
Fs = 160e6; % sample rate
R = 8;  % decimator factor
D = 1;  % differential delay
N = 5;  % number of stage
Fp = 8e6; % pass band
Fstp = 10e6; % stop band
Ap = 0.1; % attenuation in pass band
Astp = 60; % attenuation in stop band

CICDecim = dsp.CICDecimator(R, D, N);
CICCompDecim = dsp.CICCompensationDecimator(CICDecim, ...
    'DecimationFactor',2,'PassbandFrequency',Fp, ...
    'StopbandFrequency',Fstp,'SampleRate',Fs/R);
fvtool(CICDecim,CICCompDecim,...
cascade(CICDecim,CICCompDecim),'ShowReference','off','Fs',[Fs Fs/R Fs])

legend('CIC Decimator','CIC Compensator','Resulting Cascade Filter');

cicflt = cascade(CICDecim,CICCompDecim);

noise = randn(1,length(t)+1);
noise = diff(noise);
noise_cic = filter(cicflt,1,noise);

figure;plot(noise_cic)