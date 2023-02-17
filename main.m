%% 仿真单点测振仪解调过程的主函数
%%
clear;clc;close all
%% 采集设置
fs = 160e6;
ts = 1/fs;
T = 0.0015;
t = 0:ts:T-ts;
len = length(t);

fft_T = 0.001;
fft_len = fft_T*fs;
%% 测量信号
% 激光参数
fc = 40e6;              % 载波频率 Hz
lamada = 6328e-10;      % 波长
mf = 4*pi/lamada;       % 调制指数

% 振动信号参数
fv = 10e3;                 % 振动频率 Hz
vm = 0.005;               % 速度幅值 m/s
am = vm/2/pi/fv;          % 位移幅值
v = vm*cos(2*pi*fv*t);
s = am*sin(2*pi*fv*t);    % 位移信号
s = s-mean(s);
signal = sin(2*pi*fc*t+mf*s);

%% 添加宽带噪声和低频噪声
signal = awgn(signal,30,'measured');

figure;
subplot(211);plot(t,v);title('振动速度信号')
hold on;
subplot(212);plot(t,signal);title('多普勒信号')

[S_v,f_v] = Func_pufenxi(v(1:fft_len),fs,1);
[S_signal,f_signal] = Func_pufenxi(signal(1:fft_len),fs,1);
figure;
subplot(211);plot(f_v,S_v);title('振动速度信号')
hold on;
subplot(212);plot(f_signal,S_signal);title('多普勒信号')

%% 本振信号使用理想采样信号
sin_rpt = repmat([0 1 0 -1],1,length(signal)/4);
cos_rpt = repmat([1 0 -1 0],1,length(signal)/4);
[S_v,f_v] = Func_pufenxi(cos_rpt(1:fft_len),fs,1);
figure;
subplot(211);plot(f_v,S_v);title('sin')
%% 混频后滤波
ss =  signal.*sin_rpt;
sc =  signal.*cos_rpt;
 
%% cic滤波
%        (1-z^RD)^M
% H(z) = ————    *(1/RD)^M
%        (1-z^-1)^M
M = 5;     % 级联个数，改变衰减 astop = 13*M dB
R = 4000;  % 改变滤波后的采样频率 fs2 = fs/R;
D = 1;
s_sin_flt = Func_cic(ss,M,R,D);
s_cos_flt = Func_cic(sc,M,R,D);

%% 反正切
s_atan = atan2(s_sin_flt,s_cos_flt);
%% 滤波
