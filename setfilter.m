%% Fs_160M.m 仿真VD-02解调，
% 针对采集频率为160MHz
% 载波fc=40MHz
% 目标振动0-30MHz
% 速度最大值为12m/s
% 在频带范围内计算各个频点的响应
clear
clc
close all
%% 采集设置 
fs = 160e6;
ts = 1/fs;
T = 0.015;
t = 0:ts:T-ts;
len = length(t);

fpass = 35e6;
fstop = 40e6;

Num = GetLowFilterCoef(fs,fpass,fstop,0.01,100,20);
% fvtool(Num,1);
Numdelay = grpdelay(Num);
delay = floor(2*Numdelay(1))+1;

fv = 30e6;

fc = 40e6;             % 载波频率
lamada = 6328e-10;     % 波长
mf = 4*pi/lamada;      % 调制指数
vmax = 5e-3;             % 速度幅值
am = vmax/2/pi/fv;     % 位移幅值
s = am*sin(2*pi*fv*t); % 位移信号
s_sin = sin(mf*s);

[Signal_abs,f_signal2] = Func_pufenxi(s_sin,fs,0);
figure;plot(f_signal2,Signal_abs);title('sin信号频谱')

s_cos = cos(mf*s);

[Signal_abs,f_signal2] = Func_pufenxi(s_cos,fs,0);
figure;plot(f_signal2,Signal_abs);title('cos信号频谱')

s_sin_flt = filter(Num,1,s_sin);
s_cos_flt = filter(Num,1,s_cos);

phat0 = atan2(s_sin_flt,s_cos_flt);
phat1 = unwrap(phat0);
phat1 = phat1(delay:delay+0.01*fs);
vv = diff(phat1/mf)/ts;
figure
plot(vv)
title('解调信号波形')

[vv_abs,f_vv] = Func_pufenxi(vv,fs,0);
figure;plot(f_vv,vv_abs);title('解调信号频谱')

%% 目标
error = zeros(1,2501);
error_dB = zeros(1,2501);
v22 = zeros(1,2501); 

for i = 2500:2500
    
    fv = i*1000;
    
    % fv = 2.5e6;            % 振动频率
    vmax = 10;             % 速度幅值
    am = vmax/2/pi/fv;     % 位移幅值
    s = am*sin(2*pi*fv*t); % 位移信号
    v = diff(s)/ts;        % 速度信号
    %% 多普勒信号
    fc = 40e6;             % 载波频率
    lamada = 6328e-10;     % 波长
    mf = 4*pi/lamada;      % 调制指数
    signal = cos(2*pi*fc*t+mf*s);  % 已调信号
    % load('Num_band_30-34.5_45.5-50MHz_160M_0.01-80.mat')
    % signal = filter(Num,1,signal0);

    [Signal_abs,f_signal2] = Func_pufenxi(signal,fs,0);
    figure;plot(f_signal2,Signal_abs);title('已调信号频谱')

    N = fv/fs*len;  % fv的n次谐波频点对应的点数
    b = Signal_abs(1:N:end);
    Jn = floor(fc/fv);  % 0阶到Jn-1阶
    
    Firsthalf = 2:Jn;
    Secondhalf = Firsthalf+Jn;
    
    b1 = b(Firsthalf);
    b2 = flip(b(Secondhalf));
    b22 = [b1;b2];figure;subplot(211);stem(b1);hold on;stem(b2,'*')
    b3 = b2-b1;subplot(212);plot(b3);title('误差')

    %% 解调-混频
    s_sin = signal.*sin(2*pi*fc*t);

    [Ssin_abs,f_sin2] = Func_pufenxi(s_sin,fs,1);
    figure;plot(f_sin2,Ssin_abs);title('正交混频后信号频谱')

    s_cos = signal.*cos(2*pi*fc*t);

    [Scos_abs,f_cos2] = Func_pufenxi(s_cos,fs,1);
    figure;plot(f_cos2,Scos_abs);title('同相混频后信号频谱')

    %% 解调-滤波+反正切
%     滤波器系数
    if fs == 8e6
        fpass = 1.625e6;
        fstop = 2e6;
    else 
        fpass = 34.5e6;
        fstop = 40e6;
    end

    Num = GetLowFilterCoef(fs,fpass,fstop,0.01,100,20);
    % fvtool(Num,1);
    Numdelay = grpdelay(Num);
    delay = floor(Numdelay(1))+1;

    % 同向和正交混频
    s_sin_flt = filter(Num,1,s_sin);
    s_cos_flt = filter(Num,1,s_cos);

    % 解相位和解缠
    phat2 = atan2(s_sin_flt,s_cos_flt);
    phat2 = unwrap(phat2);

    % 位移和速度
    s1 = phat2/(4*pi/lamada);
    s2 = s1(delay:end);
    v2 = diff(s2)/ts;

    %% 图像
    % 输入信号
    figure;plot(s);title('输入位移信号')
    figure;plot(v);title('输入速度信号')

    % 输出信号
    figure;plot(s2);title('输出位移信号')
    
    figure;plot(v2);title('输出速度信号')
    
    [S2_abs,f_S2] = Func_pufenxi(s2,fs,0);
    figure;plot(f_S2,S2_abs);title('输出位移信号频谱')
    % 
    [V2_abs,f_V2] = Func_pufenxi(v2,fs,0);
    V2_max = max(V2_abs);
    v22(i) = V2_max;
    error(i) = vmax-V2_max;
    
    error_dB(i) = 20*log10(V2_max/vmax);
    
    figure;plot(f_V2,V2_abs);title('输出速度信号频谱')

    %% 理论计算值
    beta = mf*am;
    a = beisaier([0:0.001:200],beta,22,1);
%     
%     beta = mf*(3/2/pi/fv);
%     a = beisaier([0:0.001:200],beta,Jn);
    
end

20*log10(V2_max/vmax)

