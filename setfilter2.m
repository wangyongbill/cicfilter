%% 通过均值后直接进行反正切
% 针对采集频率为160MHz 载波fc=40MHz
% 目标振动0-30MHz 速度最大值为12m/s 在频带范围内计算各个频点的响应
clear
clc
close all
%% 采集设置
fs = 160e6;
ts = 1/fs;
T = 0.015;
t = 0:ts:T-ts;
len = length(t);
%% 噪声
fpass = 1e6;
fstop = 2e6;
low_flt = GetLowFilterCoef(fs,fpass,fstop,1,50,20);

low_noise = 0.01*filter(low_flt,1,randn(1,length(t)));
[v_abs,f_v] = Func_pufenxi(low_noise,fs,0);
figure;plot(f_v,v_abs);title('低频噪声频谱')

%%
fv = 10e3;             % 振动频率
fc = 40e6;             % 载波频率
lamada = 6328e-10;     % 波长
mf = 4*pi/lamada;      % 调制指数
vm = 1;                % 速度幅值

v = vm*cos(2*pi*fv*t); % 速度信号
figure;plot(t,v);title("振动速度波形")

[v_abs,f_v] = Func_pufenxi(v,fs,0);
figure;plot(f_v,v_abs);title('振动速度信号频谱')

s = cumsum(v)/fs;
s = s - mean(s);
figure;plot(t,s);title("振动位移波形")

[s_abs,f_s] = Func_pufenxi(s(1:0.001*fs),fs,0);
figure;plot(f_s,s_abs);title('振动位移信号频谱')
%% 调制到40M的载波上
signal = cos(2*pi*fc*t+mf*s)+low_noise;   % 已调信号 低频噪声0-1M带宽
% figure;plot(t,signal);title("带噪调制信号波形")
%
% [signal_abs,f_signal] = Func_pufenxi(signal,fs,0);
% figure;plot(f_signal,signal_abs);title('带噪调制信号频谱')
%% 带通滤波
% signal = filter(band_filter,1,signal);
signal = awgn(signal,30,'measured');
figure;plot(t,signal);title("滤波后调制信号波形")
[signal_abs,f_signal] = Func_pufenxi(signal,fs,0);
figure;plot(f_signal,signal_abs);title('滤波后调制信号频谱')

%% 解调

s_sin = signal.*repmat([0 1 0 -1],1,length(signal)/4);
[Ssin_abs,f_sin2] = Func_pufenxi(s_sin,fs,0);
figure;plot(f_sin2,Ssin_abs);title('正交混频后信号频谱')

s_cos = signal.*repmat([1 0 -1 0],1,length(signal)/4);
[Scos_abs,f_cos2] = Func_pufenxi(s_cos,fs,0);
figure;plot(f_cos2,Scos_abs);title('同相混频后信号频谱')

% 同向和正交混频
sumNum = 16;
noise_Mat = reshape(s_sin,[sumNum,length(s_sin)/sumNum]);
s_sin_flt = sum(noise_Mat,1)/sumNum;
[Ssin_flt_abs,f_flt_sin2] = Func_pufenxi(s_sin_flt,fs/sumNum,0);
figure;plot(f_flt_sin2,Ssin_flt_abs);title('正交混频滤波后信号频谱')

noise_Mat = reshape(s_cos,[sumNum,length(s_cos)/sumNum]);
s_cos_flt = sum(noise_Mat,1)/sumNum;
[Scos_flt_abs,f_flt_cos2] = Func_pufenxi(s_cos_flt,fs/sumNum,0);
figure;plot(f_flt_cos2,Scos_flt_abs);title('同相混频滤波后信号频谱')

% 解相位和解缠
phat2 = atan2(s_sin_flt,s_cos_flt);
phat2 = unwrap(phat2);

% 位移和速度
s1 = phat2/(4*pi/lamada);

v20 = diff(s1)/ts/sumNum;
v2 = v20(end-0.001*fs/sumNum+1:end);

%%
if 1
    % 输出信号
    figure;plot(v2);title('输出速度信号')
    [V2_abs,f_V2] = Func_pufenxi(v2,fs/sumNum,0);
    figure;plot(f_V2,V2_abs);title('输出速度信号频谱')
    
    s2 = cumsum(v2)/fs*sumNum;
    figure;plot(s2);title('输出位移信号')
    [S2_abs,f_S2] = Func_pufenxi(s2,fs/sumNum,0);
    figure;plot(f_S2,S2_abs);title('输出位移信号频谱')
end

if 0
    fpass = 4e6;
    fstop = 8e6;
    low_flt = GetLowFilterCoef(fs,fpass,fstop,0.01,100,20);
    
    v2 = filter(low_flt,1,v2);
    
    lag = 8;
    % 输出信号
    figure;plot(v2(1:lag:end));title('输出速度信号')
    [V2_abs,f_V2] = Func_pufenxi(v2(1:lag:end),fs/lag,0);
    figure;plot(f_V2,V2_abs);title('输出速度信号频谱')
    
    s2 = cumsum(v2(1:lag:end))/fs*lag;
    figure;plot(s2(1:end));title('输出位移信号')
    [S2_abs,f_S2] = Func_pufenxi(s2(1:end),fs/lag,0);
    figure;plot(f_S2,S2_abs);title('输出位移信号频谱')
end