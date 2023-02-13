function [S,f] = Func_pufenxi(y,fs,halfscale)
% ���룺
%    y�������ź�
%    fs:����Ƶ��
%    halfscale: 0,�򻭳�0-fs��Ƶ�Σ�1���򻭳�0-fs/2��Ƶ��
% �ó�������ݽ����׷�����������ݵ�Ƶ�ף������굥λΪHz
% ������Ϊ���������һ����(0Hz)���ӱ�����2����(N+1)/2����ӱ�
% ������Ϊż�������һ����(OHz)���ӱ�����2����N/2����ӱ�����N/2+1��(N/2Hz)���ӱ�
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



