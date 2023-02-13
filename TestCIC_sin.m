%% 设计降采样滤波器，测试滤波器的效果
% 测试CIC有效性
clear
clc
close all
%% 采集设置
fs = 160e6;
ts = 1/fs;
T = 0.00015;
t = 0:ts:T-ts;
len = length(t);
%% 信号生成
fv1 = 1e4;             % 振动频率
fv2 = 0;             % 振动频率
fv3 = 0;             % 振动频率
fv4 = 0;             % 振动频率
vm = 0.01;                % 速度幅值 m/s

v = vm*cos(2*pi*fv1*t)+vm*cos(2*pi*fv2*t)+vm*cos(2*pi*fv3*t)+vm*cos(2*pi*fv4*t); % 速度信号
figure;plot(t,v);title("振动速度波形")
[v_abs,f_v] = Func_pufenxi(v,fs,0);
figure;plot(f_v,v_abs);title('振动速度信号频谱')

%% 测量信号
fc = 40e6;
lamada = 6328e-10;
mf = 4*pi/lamada;
s = cumsum(v)/fs;
s = s-mean(s);
signal = cos(2*pi*fc*t+mf*s);
signal = awgn(signal,30,'measured');
figure;plot(t,signal);title("带噪信号波形")
[v_abs,f_v] = Func_pufenxi(signal,fs,0);
figure;plot(f_v,v_abs);title('带噪信号频谱')
%% 混频
s_sin = signal.*repmat([0 1 0 -1],1,length(signal)/4);
s_cos = signal.*repmat([1 0 -1 0],1,length(signal)/4);
%% 解调下采样频率，cic下采样
R = 4000;
D = 1;
M = 1;

Hd = dsp.CICDecimator( ...
    'DecimationFactor', R, ...
    'DifferentialDelay', D, ...
    'NumSections', M);

numPoint = length(signal)/10*R;
src1 = dsp.SignalSource(s_sin',numPoint);
src2 = dsp.SignalSource(s_cos',numPoint);

for i = 1:2048*R
    s_sin_cic((i-1)*numPoint/R+1:i*numPoint/R) = hcic(src1());
end
for i = 1:2048*R
    s_cos_cic((i-1)*numPoint/R+1:i*numPoint/R) = hcic(src2());
end

[CIC_abs,f_CIC] = Func_pufenxi(s_sin_cic,fs/R,0);
figure;plot(f_CIC,CIC_abs);title('cic滤波后频谱')

[CIC_abs,f_CIC] = Func_pufenxi(s_cos_cic,fs/R,0);
figure;plot(f_CIC,CIC_abs);title('cic滤波后频谱')
%%
phat2 = atan2(s_sin_cic,s_cos_cic);
phat2 = unwrap(phat2);

% 位移和速度
s1 = phat2/(4*pi/lamada);
v2 = [0 diff(s1)*fs/R];

%% cic
v3 = v2(end-fs/R*0.0001+1:end);
figure;plot(v3);title('输出速度信号')
[V3_abs,f_V3] = Func_pufenxi(v3,fs/R,0);
figure;plot(f_V3,V3_abs);title('输出速度信号频谱')

%% 
