clear
clc
close all

fs = 160e6;
ts = 1/fs;
T = 0.0015;
t = 0:ts:T-ts;
len = length(t);

x = sin(2*pi*10e3*t);
figure;plot(x)

figure;plot(y)

R = 1000;
D = 1;
b = [1 zeros(1,R*D-1) -1];    % 抽取倍数等于差分时延
a = [1 -1];
Hd = Func_designCIC(R,D,1);
fvtool(Hd)

function Hd = Func_designCIC(R,D,M)

decf    = R;     % Decimation Factor
diffd   = D;     % Differential Delay
numsecs = M;     % Number of Sections

Hd = dsp.CICDecimator( ...
    'DecimationFactor', decf, ...
    'DifferentialDelay', diffd, ...
    'NumSections', numsecs);
end
% [EOF]
