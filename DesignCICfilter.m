clear, clc, close all
% 设置补偿滤波器的阶数
FIR_Order = 220; % 滤波器阶数比滤波器系数个数少1
% 频谱的频率点数
Num_Freq_Points = 512*512; % Also try, 256 or 1024
% 设置显示的最小幅值(dB)
Threshold = -200; % dB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define parameters of your desired CIC filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 5;     % 滤波器级联级数
R = 4000;  % 抽取因子
N = 1;     % 差分时延
Fs = 160e6;    % 原始采样频率
Fp = 4200;     % 通带截止频率(-6dB)

% 检查通带范围，需不大于降采样后采样频率的一半
if Fp >= 0.5*Fs/R 
    disp(['Fp必须小于Fs/R = ',num2str(Fs/R) ,' Hz的1/2'])
    return 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CIC滤波器的幅频
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 归一化频率范围
Freq = 0: 1/(2*Num_Freq_Points): 0.5-1/(2*Num_Freq_Points);
% CIC幅频函数
Spec_CIC = (R*N*sin(pi*Freq*R*N)./sin(pi*Freq)).^M;
Spec_CIC(1) = (R*N)^M;
Spec_Mag_CIC = abs(Spec_CIC); % 幅度
Spec_Mag_CIC = Spec_Mag_CIC/max(Spec_Mag_CIC); % 归一化幅度
Spec_Mag_dB_CIC = 20*log10(Spec_Mag_CIC); % 转换到dB
figure(1), clf
Freq_Plot_Axis_1 = (0:Num_Freq_Points-1)*Fs/(2*Num_Freq_Points);
plot(Freq_Plot_Axis_1, Spec_Mag_dB_CIC, '-k')
% axis([0, Fs/2, Threshold, 5])
xlim([0, Fs/2])
xlabel('Hz'), ylabel('dB'),
title('CIC 幅度响应'), grid on, zoom on
% Determine & plot CIC filter aliasing after decimation
Mainlobe_Indices = 1:round(2*Num_Freq_Points/(R));
Mainlobe = Spec_Mag_dB_CIC(Mainlobe_Indices);
Mainlobe_Flipped = fliplr(Mainlobe);
% Find "aliased-mainlobe < mainlobe" freq indices
Indices = [1]; % Intialize
for Loop = 2:round(2*Num_Freq_Points/(R))
    if Mainlobe(Loop) >= Mainlobe_Flipped(Loop)
        Indices = [Indices, Loop];
    else
        end
end
hold on
plot(Freq_Plot_Axis_1(Indices), Mainlobe_Flipped(Indices), '-r')
hold off
text(Fs/R, 0, 'Red shows the most significant')
text(Fs/R, -7, ' (but not the total) CIC filter')
text(Fs/R, -14, ' aliasing after decimation')

% End of CIC filter design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cic的fir补偿滤波器设计

Num_FIR_Design_Freq_Points = 4*FIR_Order;
FIR_Freq_Vector = linspace(0, 1, Num_FIR_Design_Freq_Points);

% 生成cic幅频幅值函数的倒数，作为补偿滤波器的幅度
CIC_Freq_vector = 0:round(Num_FIR_Design_Freq_Points*Fp/(Fs/(2*R)));
CIC_Freq_vector = (1/R)*CIC_Freq_vector/Num_FIR_Design_Freq_Points;
% 根据公式，生成倒数
Inverse_CIC = 1./((1/(R*N)^M)*(sin(pi*CIC_Freq_vector*R*N/2)...
    ./sin(pi*CIC_Freq_vector/2)).^M);
%  Richard Lyons – March, 2020 Page 16
Inverse_CIC(1) = 1; % Eliminate the NaN first sample
Inverse_CIC(Num_FIR_Design_Freq_Points) = 0;
Freq_Plot_Axis_2 = (0:Num_FIR_Design_Freq_Points-1)*(Fs/(2*R))...
    /(Num_FIR_Design_Freq_Points);
figure(2), clf
plot(Freq_Plot_Axis_2, Inverse_CIC, '-bs', 'markersize', 2)
title('Desired FIR resp.= Blue, Actual FIR resp.= Red')
ylabel('Linear'), xlabel('Hz'), grid on, zoom on
% Compute FIR compensation filter coeffs
FIR_Coeffs = fir2(FIR_Order, FIR_Freq_Vector, ...
    Inverse_CIC, chebwin(FIR_Order+1, 50));
disp(['FIR coefficients = ', num2str(FIR_Coeffs)])
% Compute freq response of a unity-gain FIR Comp. filter
Num_Freq_Points = 512;
Spec_Comp_Filter = fft(FIR_Coeffs, 2*Num_Freq_Points);
Spec_Comp_Filter = Spec_Comp_Filter(1:Num_Freq_Points); % Pos. freqs only
Mag_Comp_Filter = abs(Spec_Comp_Filter);
Mag_dB_Comp_Filter = 20*log10(Mag_Comp_Filter);
%Freq_Plot_Axis_3 = (0:Num_Freq_Points-1)*(Fs/R)/(Num_Freq_Points);
Freq_Plot_Axis_3 = (0:Num_Freq_Points-1)*(Fs/(2*R))/(Num_Freq_Points);
figure(2)
hold on
plot(Freq_Plot_Axis_3, Mag_Comp_Filter, '-r', 'markersize', 2)
hold off
axis([0, Fs/(2*R), 0, 1.2*max(Inverse_CIC)])
figure(3), clf
plot(Freq_Plot_Axis_3, Mag_dB_Comp_Filter)
axis([0, Fs/(2*R), -80, max(Mag_dB_Comp_Filter)+5])
xlabel('Hz'), ylabel('dB')
title('FIR Comp. Filter mag. resp.'), grid on, zoom on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of FIR compensation filter design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiply the unity-gain FIR comp. filter complex freq response times
% the CIC filter's complex **mainlobe-only** freq response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute CIC filter's mainlobe-only response (using Num_Freq_Points)
CIC_Mainlobe_Freq_vector = 0:Num_Freq_Points -1;
CIC_Mainlobe_Freq_vector = 0.5*(2/R)*...
    CIC_Mainlobe_Freq_vector/Num_Freq_Points;
CIC_Mainlobe = ((1/(R*N))*(sin(pi*CIC_Mainlobe_Freq_vector*R*N/2)...
    ./sin(pi*CIC_Mainlobe_Freq_vector/2))).^M;
CIC_Mainlobe(1) = 1; % Zero Hz

% Compute the frequency magnitude response of the cascaded filters
Cascaded_Filter_Spec_Mag = abs(CIC_Mainlobe...
    .*Spec_Comp_Filter);
%  Richard Lyons – March, 2020 Page 17

% Compute cascaded response in dB for plotting
Cascaded_Filter_Spec_Mag_dB = 20*log10(Cascaded_Filter_Spec_Mag);
Freq_Plot_Axis_3 = (0:Num_Freq_Points-1)*(Fs/(2*R))/(Num_Freq_Points);
% Plot cascaded filters' combined freq. magnitude response
figure(4), clf
subplot(2,1,1)
plot(Freq_Plot_Axis_3, Cascaded_Filter_Spec_Mag_dB, '-k')
hold on
plot(Freq_Plot_Axis_1(1:51), Spec_Mag_dB_CIC(1:51), '-r')
plot(Freq_Plot_Axis_3, Mag_dB_Comp_Filter, '-b')
hold off
axis([0, Fs/(2*R), Threshold, max(Mag_dB_Comp_Filter)+10])
xlabel('Hz'), ylabel('dB')
title('Black = cascaded filter, Blue = comp. filter, Red = CIC')
grid on, zoom on

subplot(2,1,2)
plot(Freq_Plot_Axis_3, Cascaded_Filter_Spec_Mag_dB, '-k')
hold on
plot(Freq_Plot_Axis_1(1:51), Spec_Mag_dB_CIC(1:51), '-r')
plot(Freq_Plot_Axis_3, Mag_dB_Comp_Filter, '-b')
hold off
title('Black = cascaded filter, Blue = comp. filter, Red = CIC')
axis([0, 1.2*Fp, -6-N, max(Mag_dB_Comp_Filter)+2])
xlabel('Hz'), ylabel('dB')
grid on, zoom on
