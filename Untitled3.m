close all
clear all
clc

Fs = 160e6; % sample rate
R = 20;  % decimator factor
D = 1;  % differential delay
Fp = 2e6; % pass band
Fstp = 4e6; % stop band
Ap = 0.1; % attenuation in pass band
Astp = 60; % attenuation in stop band
hcic = design(fdesign.decimator(R,'cic',D, Fp, Astp, Fs),'SystemObject',true);
cic_comp = design(fdesign.ciccomp(hcic.DifferentialDelay, ...
            hcic.NumSections,Fp,Fstp,Ap,Astp,Fs/R), 'SystemObject',true);
fvtool(hcic,cic_comp,...
cascade(hcic,cic_comp),'ShowReference','off','Fs',[Fs Fs/R Fs])
legend('CIC Decimator','CIC Compensator','Resulting Cascade Filter');

%  
% %set system parameter
f = 30e3;    %The frequency of the local oscillator signal
% Fs = 160e6;   %sampling frequency
% N =  24;         %Quantitative bits
% L = 1600000;
%  
% %Generating an input signal
% t =0:1/Fs:(1/Fs)*(L-1);          %Generating the time series of sampling frequencies
% sc =sin(2*pi*500e3*t)+sin(2*pi*f*t);        %a sinusoidal input signal that produces a random starting phase
%  
% b =[1,-1];%integerator
% a =[1,-1];%comb 
%  
% %integerator
% c1=filter(1,b,sc);
% c2=filter(1,b,c1);
% c3=filter(1,b,c2);
%  
% % y = downsample(c3,10);
% y = c3(1:1000:end);
% %comb
% i1 =filter(a,1,y);
% i2 =filter(a,1,i1);
% i3 =filter(a,1,i2);
% sf = i3;
% % sf = i3(1:1000:end);
% % sf = downsample(i3,1000);
% 
% figure(1),
% plot(sf);
% title('sf','fontsize',8);
% [S,f] = Func_pufenxi(sf,Fs/1000,0);
% figure;plot(f,S/max(S))
%  
% %% 等效3个fir级联后抽取
% % N = 3 R = 4 M = 1
% coe = [1 1 1 1];
% coe1 = conv(coe,coe);
% coe2 = conv(coe1,coe);
% sf1 =filter(coe2,1,sc);
% sf2 = downsample(sf1,4);
%  
% %% 
% figure(1),
% subplot(211),stem(t(1:1600),sc(1:1600));
% xlabel('时间(t)','fontsize',8);
% ylabel('幅度(dB)','fontsize',8);
% title('sc','fontsize',8);
% subplot(212),stem(t(1:400),sf(1:400));
% xlabel('时间/4(t)','fontsize',8);
% ylabel('幅度(dB)','fontsize',8);
% title('sf','fontsize',8);
% 
% figure(2),
% subplot(211),stem(t(1:1600),sc(1:1600));
% xlabel('时间(t)','fontsize',8);
% ylabel('幅度(dB)','fontsize',8);
% title('sc','fontsize',8);
% subplot(212),stem(t(1:400),sf2(1:400));
% xlabel('时间/4(t)','fontsize',8);
% ylabel('幅度(dB)','fontsize',8);
% title('sf','fontsize',8);
% 
% figure(3)
% plot(sf2 - sf) 
% 
% 
% % % Filename: CIC_compensation_filter.m
clear
clc
% 补偿滤波器阶数，阶数加一等于滤波器系数个数
FIR_Order = 31; 
% 谱线数
Num_Freq_Points = 512; %  256 or 1024
% Set minimum mag. level (in dB) in spectral mag. plots
Threshold = -80; % Minimum dB plot level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameters of your desired CIC filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 1; % Number of cascaded CIC stages.
R = 4; % Decimation sample rate-change factor.
N = 1; % Differential delay after sample rate change.
Fs = 1000; % CIC input sample rate before decimation (Hz).
Fp = 40; % Final cascaded-filters' passband cutoff freq
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
[Taps, W, beta, ftype] = firpmord([Fp, 1.3*Fp], [1,0], [0.05, 0.01], Fs);
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
