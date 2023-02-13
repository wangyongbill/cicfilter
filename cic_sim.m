close all
clc;
clear;
%% 正弦叠加信号
sin_f0 = 30e3;
sin_f1 = 500e3;
sin_f3 = 60e6;

Fs = 160e6;
N = 16e4;
n = 0:N-1;
time = n/Fs;

y0 = sin(2*pi*sin_f0.*time);
y1 = sin(2*pi*sin_f1.*time);
y3 = sin(2*pi*sin_f3.*time);
y = y0+y1+y3;

base_len = length(y);
f_label = (0:base_len-1)*Fs/base_len;

downNum = 1000;
delayNum = 8;

y_len = base_len/downNum;
y_label = (0:y_len-1)*Fs/downNum/y_len;

base_sig = y;
base_fft = abs(fft(base_sig));
figure(1);subplot(611);plot(base_sig(1:1024),'b');grid;title('原始时域信号')

subplot(612);plot(f_label(1:end/2),base_fft(1:end/2),'b');grid;title('原始频域波形')

%% 自带函数设计的cic
CFs = 160e6;
CAstp = 0.01;
CFp = 50e3;
% cic设计
hcic = design(fdesign.decimator(downNum,'cic',delayNum,CFp,CAstp,CFs),'SystemObject',true);

NumPerFrame = 2000;  % 每帧多少个点
src = dsp.SignalSource(base_sig',NumPerFrame);
fir_out = zeros(1,N/downNum);

for i = 1:N/NumPerFrame
    fir_out((i-1)*NumPerFrame/downNum+1:i*NumPerFrame/downNum) = hcic(src());
end

fir_offt = abs(fft(fir_out));
subplot(613)
plot(y_label(1:end/2),fir_offt(1:end/2),'b');grid
title('过函数CIC滤波器频域信号')

subplot(614);
plot(fir_out(1:128),'b');grid;
title('过函数CIC滤波器时域信号')

%% 自实现的cic
cic_i1 = zeros(1,N);
cic_i1(1) = base_sig(1);
% 积分
for i = 2:N
    cic_i1(i) = cic_i1(i-1)+base_sig(i);
end
cic_i2 = cic_i1;

%% 选择等效形式,设置抽取的位置在积分后还是在差分后
%  积分后进行抽取，能降低寄存器数量
slct = 1;
if slct == 1 %% 先抽取后差分
    cic_d = zeros(1,N/downNum);
    cic_c2 = zeros(1,N/downNum);
    % 抽取
    for i = 1:N/downNum
        cic_d(i) = cic_i2((i-1)*downNum+1);
    end
    cic_c2(1) = cic_d(1);
    % 差分，差分时延为delayNum
    for i = delayNum+1:N/downNum
        cic_c2(i) = cic_d(i)-cic_d(i-delayNum);
    end
    cic_fft = abs(fft(cic_c2)/downNum/delayNum);
else  %先差分后抽取
    cic_d = zeros(1,N/downNum);
    cic_c2 = zeros(1,N);
    % 差分，差分时延为delayNum*R
    for i = delayNum*downNum+1:N
        cic_c2(i) = cic_i2(i)-cic_i2(i-delayNum*downNum);
    end
    %抽取
    for i = 1:N/downNum
        cic_d(i) = cic_i2((i-1)*downNum+1);
    end
    cic_fft = abs(fft(cic_d)/downNum/delayNum);
end
%%
subplot(615);
plot(y_label(1:end/2),cic_fft(1:end/2),'b');grid;
title('自实现CIC滤波器频域')

subplot(616)
plot(cic_d,'b');grid
title('自实现CIC滤波器时域信号')

%% cic滤波器的频域响应
figure(2)
% 积分器
b1 = [0,1];
a1 = [-1,1];
h1 = freqz(b1,a1,512);
subplot(611);plot(abs(h1),'b');grid;
title('积分器幅频响应')
subplot(612);plot(angle(h1),'b');grid;
title('积分器相频响应')

% 微分器
b2 = [-1,0,0,0,0,0,0,0,1];
a2 = [0,0,0,0,0,0,0,0,1];
h2 = freqz(b2,a2,512);
subplot(613);plot(abs(h2),'b');grid;
title('微分器幅频响应')
subplot(614);plot(angle(h2),'b');grid;
title('微分器项相频响应')

% cic
b3 = [-1,0,0,0,0,0,0,0,1];
a3 = [0,0,0,0,0,0,0,0,-1,1];
h3 = freqz(b3,a3,512);
subplot(615);plot(abs(h3),'b');grid
title('cic幅频响应')
subplot(616);plot(angle(h3),'b');grid
title('cic相频响应')
