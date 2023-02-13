function b = GetBandFilterCoef(Fs,Fstop1,Fpass1,Fpass2,Fstop2,DstopdB1,DpassdB,DstopdB2)

dens   = 20;              % Density Factor

Dpass = (10^(DpassdB/20)-1)/(10^(DpassdB/20)+1);%0.057501127785;  % Passband Ripple
if DstopdB1>0
    Dstop1 = 10^(DstopdB1/(-20));%1e-05;           % Stopband Attenuation
else
    Dstop1 = 10^(DstopdB1/(20));
end

if DstopdB2>0
    Dstop2 = 10^(DstopdB2/(-20));%1e-05;           % Stopband Attenuation
else
    Dstop2 = 10^(DstopdB2/(20));
end

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 ...
                          0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
end
