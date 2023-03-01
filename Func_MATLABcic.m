% function out = Func_Matlabcic(signal,R,D,Fp,Fstp.Ap,Astp)

Fs = 1;    % sample rate
R = 4;     % decimator factor
D = 1;     % differential delay
Fp = 0.005; % pass band
Fstp = 0.0075; % stop band
Ap = 0.1;  % attenuation in pass band
Astp = 60; % attenuation in stop band

hcic = design(fdesign.decimator(R,'cic',D, Fp, Astp, Fs),'SystemObject',true);
cic_comp = design(fdesign.ciccomp(hcic.DifferentialDelay, ...
    hcic.NumSections,Fp,Fstp,Ap,Astp,Fs/R), 'SystemObject',true);

F0 = Fs/40;
N = 1024;
t = (0:N-1)';  % length: 1024
x = fi(sin(2*pi*F0/Fs*t),true,16,15);
SamplesPerFrame = 64;
src = dsp.SignalSource(x,SamplesPerFrame);
y = zeros(N/SamplesPerFrame,SamplesPerFrame/R);
for ii = 1:length(x)/SamplesPerFrame
     y(ii,:) = hcic(src());    
end
yy = y.';
rr = yy(:);
figure;plot(yy(:))

