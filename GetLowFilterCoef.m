function [Num] = GetLowFilterCoef(Fs,Fpass,Fstop,DpassdB,DstopdB,dens)

%% 生成低通滤波器系数 
% Fs 采样频率； Fpass 通带截止频率；Fstop阻带截止频率；
% DpassdB 通带纹波；DstopdB 阻带抑制比；dens 频率密度因子

    Dpass = (10^(DpassdB/20)-1)/(10^(DpassdB/20)+1);%0.057501127785;  % Passband Ripple
    if DstopdB>0
        Dstop = 10^(DstopdB/(-20));%1e-05;           % Stopband Attenuation
    else
        Dstop = 10^(DstopdB/(20));
    end
    % Calculate the order from the parameters using FIRPMORD.
    [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

    % Calculate the coefficients using the FIRPM function.
    Num  = firpm(N, Fo, Ao, W, {dens});
end