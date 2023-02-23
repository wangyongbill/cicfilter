close all
clear
clc

% 补偿滤波器阶数，阶数加一等于滤波器系数个数
FIR_Order = 50; 
% 谱线数
Num_Freq_Points = 512; %  256 or 1024
% Set minimum mag. level (in dB) in spectral mag. plots
Threshold = -200; % Minimum dB plot level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters of your desired CIC filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 1; % Number of cascaded CIC stages.
R = 4000; % Decimation sample rate-change factor.
N = 1; % Differential delay after sample rate change.
Fs = 160e6; % CIC input sample rate before decimation (Hz).
Fp = 7.5e3; % Final cascaded-filters' passband cutoff freq
Fstop = 15e3;
% (-6 dB mag. point) after decimation (Hz).
disp(' ')
disp(['CIC filter order = ', num2str(M), ' stages'])
disp(['CIC filter input sample rate = ', num2str(Fs), ' Hz'])
disp(['CIC decimation factor R = ', num2str(R)])
disp(['CIC filter output sample rate = ', num2str(Fs/R), ' Hz'])
disp(['Final cascaded-filter -6 dB passband width = ', num2str(Fp), 'Hz'])
disp(['FIR filter number of taps = ', num2str(FIR_Order+1), ' taps'])
% Just for fun, estimate number of taps for a traditional
% tapped-delay line (non-CIC) "Parks-McClellan-designed"
% lowpass FIR filter.
[Taps, W, beta, ftype] = firpmord([Fp, Fstop], [1,0], [0.01, 0.0001], Fs);
disp(['[Equivalent tapped-delay line FIR lowpass filter' ...
    ' requires approximately ', num2str(Taps), ' taps]']), disp(' ')
% Check to see that passband width, Fp, is not too large
if Fp >= 0.5*Fs/R % Half the final output sample rate
    beep, pause(0.5), beep, disp(' ')
    disp('WARNING !'), disp('WARNING !')
    disp(['Fp = ',num2str(Fp),' Hz passband width is too large!'])
    disp(['Fp must be less than one half the final Fs/R = ',...
        num2str(Fs/R) ,' Hz sample rate!'])
    %  Richard Lyons – March, 2020 Page 15
    return % Stop all processing
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot target CIC filter positive-freq. mag. response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define frequency vector (Freq) in the range of 0 -to- 0.5
Freq = 0:1/(2*Num_Freq_Points):0.5 -1/(2*Num_Freq_Points);
% CIC filter frequency response equation
Spec_CIC = (sin(pi*Freq*R*N)./sin(pi*Freq)).^M;

Spec_CIC(1) = 1; % Set correct DC (zero Hz) gain
Spec_Mag_CIC = abs(Spec_CIC); %Magnitude
Spec_Mag_CIC = Spec_Mag_CIC/max(Spec_Mag_CIC); % Set max mag. to one
Spec_Mag_dB_CIC = 20*log10(Spec_Mag_CIC); % Decibels
figure(1), clf
Freq_Plot_Axis_1 = (0:Num_Freq_Points-1)*Fs/(2*Num_Freq_Points);
plot(Freq_Plot_Axis_1, Spec_Mag_dB_CIC, '-k')
axis([0, Fs/2, Threshold, 5])
xlabel('Hz'), ylabel('dB'),
title('CIC mag resp. before decimation'), grid on, zoom on
% Determine & plot CIC filter aliasing after decimation
Mainlobe_Indices = 1:round(2*Num_Freq_Points/(R));
Mainlobe = Spec_Mag_dB_CIC(Mainlobe_Indices);
Mainlobe_Flipped = fliplr(Mainlobe);
% Find "aliased-mainlobe < mainlobe" freq indices
Indices = [1]; % Intialize
for Loop = 2:round(2*Num_Freq_Points/(R))
    if Mainlobe(Loop) >= Mainlobe_Flipped(Loop)
        Indices = [Indices, Loop];
    else, end
end
hold on
plot(Freq_Plot_Axis_1(Indices), Mainlobe_Flipped(Indices), '-r')
hold off
text(Fs/R, 0, 'Red shows the most significant')
text(Fs/R, -7, ' (but not the total) CIC filter')
text(Fs/R, -14, ' aliasing after decimation')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of CIC filter design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design the FIR compensation filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Num_FIR_Design_Freq_Points = 4*FIR_Order;
FIR_Freq_Vector = linspace(0, 1, Num_FIR_Design_Freq_Points);
% Create freq vector for defining inverse CIC mag. response
% over the FIR compensation filter's 'Fp' passband.
CIC_Freq_vector = 0:round(Num_FIR_Design_Freq_Points*Fp/(Fs/(2*R)));
CIC_Freq_vector = (1/R)*CIC_Freq_vector/Num_FIR_Design_Freq_Points;
% Define an inverse CIC magnitude vector
Inverse_CIC = 1./((1/(R*N)^M)*(sin(pi*CIC_Freq_vector*R*N/2)...
    ./sin(pi*CIC_Freq_vector/2)).^M);
% Inverse_CIC = 1./((1/(R*N)^M)*(sin(pi*CIC_Freq_vector*R*N/2)...
%     ./sin(pi*CIC_Freq_vector/2)).^M);
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
