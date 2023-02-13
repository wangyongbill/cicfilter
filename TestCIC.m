%% 设计降采样滤波器，测试滤波器的效果
% 针对采集频率为160MHz 载波fc=40MHz
% 对比：解调低通

clear
clc
close all
%% 采集设置
fs = 160e6;
ts = 1/fs;
T = 0.015;
t = 0:ts:T-ts;
len = length(t);
%% 振动信号生成
fv = 1e3;             % 振动频率
fc = 40e6;             % 载波频率
lamada = 6328e-10;     % 波长
mf = 4*pi/lamada;      % 调制指数
vm = 1;                % 速度幅值 m/s

v = vm*cos(2*pi*fv*t); % 速度信号
figure;plot(t,v);title("振动速度波形")

[v_abs,f_v] = Func_pufenxi(v,fs,0);
figure;plot(f_v,v_abs);title('振动速度信号频谱')

s = cumsum(v)/fs;
s = s-mean(s);
figure;plot(t,-s);title("振动位移波形")

[s_abs,f_s] = Func_pufenxi(s(1:0.001*fs),fs,0);
figure;plot(f_s,s_abs);title('振动位移信号频谱')

%% 调制到40M的载波上
signal = cos(2*pi*fc*t+mf*s);   % 已调信号
figure;plot(t,signal);title("带噪调制信号波形")

[signal_abs,f_signal] = Func_pufenxi(signal,fs,0);
figure;plot(f_signal,signal_abs);title('带噪调制信号频谱')
%% 混频
SS = repmat([0 1 0 -1],1,length(signal)/4);
s_sin = signal.*SS;

[Ssin_abs,f_sin2] = Func_pufenxi(s_sin,fs,0);
figure;plot(f_sin2,Ssin_abs);title('正交混频后信号频谱')

SC = repmat([1 0 -1 0],1,length(signal)/4);
s_cos = signal.*SC;

[Scos_abs,f_cos2] = Func_pufenxi(s_cos,fs,0);
figure;plot(f_cos2,Scos_abs);title('同相混频后信号频谱')

%% 下采样
s_sin_flt = filter(Num,1,s_sin);
[Ssin_flt_abs,f_flt_sin2] = Func_pufenxi(s_sin_flt,fs,0);
figure;plot(f_flt_sin2,Ssin_flt_abs);title('正交混频滤波后信号频谱')

s_cos_flt = filter(Num,1,s_cos);
[Scos_flt_abs,f_flt_cos2] = Func_pufenxi(s_cos_flt,fs,0);
figure;plot(f_flt_cos2,Scos_flt_abs);title('同相混频滤波后信号频谱')

% 解相位和解缠
phat2 = atan2(s_sin_flt,s_cos_flt);
phat2 = unwrap(phat2);

% 位移和速度
s1 = phat2/(4*pi/lamada);

v20 = diff(s1)/ts;
v2 = v20(end-0.001*fs+1:end);
lag = 2;
