function [S,f] = Func_pufenxi(y,fs,halfscale)
% 输入：
%    y：输入信号
%    fs:采样频率
%    halfscale: 0,则画出0-fs的频段；1，则画出0-fs/2的频段
% 该程序对数据进行谱分析，输出数据的频谱，横坐标单位为Hz
% 若点数为奇数，则第一个点(0Hz)不加倍，第2到第(N+1)/2个点加倍
% 若点数为偶数，则第一个点(OHz)不加倍，第2到第N/2个点加倍，第N/2+1点(N/2Hz)不加倍
% 
length_y = length(y);

Y = fft(y);
Y_abs = abs(Y);
Y_Am = Y_abs;

if halfscale == 0    
    
    Y_Am(1) = Y_Am(1)/length_y;
    disfreqlen = length_y;   
    
    if (mod(length_y,2) == 0)                 
        Y_Am(2:disfreqlen/2-1) = 2*Y_Am(2:disfreqlen/2-1)/length_y;
        Y_Am(disfreqlen/2) = Y_Am(disfreqlen/2)/length_y;   
        Y_Am(disfreqlen/2+1:end) = 2*Y_Am(disfreqlen/2+1:end)/length_y;
    else    
        Y_Am(2:disfreqlen) = 2*Y_Am(2:disfreqlen)/length_y;
    end   
    
    S = Y_Am;
    f = (0:length_y-1)/length_y*fs;
else 
    if (mod(length_y,2) == 0)   
        disfreqlen = length_y/2+1;
        Y_Am(1) = Y_Am(1)/length_y;
        Y_Am(2:disfreqlen-1) = 2*Y_Am(2:disfreqlen-1)/length_y;
        Y_Am(disfreqlen) = Y_Am(disfreqlen)/length_y;   
    else 
        disfreqlen = (length_y+1)/2;
        Y_Am(1) = Y_Am(1)/length_y;
        Y_Am(2:disfreqlen) = 2*Y_Am(2:disfreqlen)/length_y;
    end   
    S = Y_Am(1:disfreqlen);
    f = (0:disfreqlen-1)/length_y*fs;
end



