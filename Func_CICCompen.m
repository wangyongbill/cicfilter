function cic_comp_flt1 = Func_CICCompen(Fpass,Fstop,Aps,Asp,K,Rate,Delay,fs)
%% CIC滤波器
D = Rate;      %% 抽取因子
M = Delay;     %% 抽取后的差分时延   主瓣带宽：2*pi/D*M
N = D*M;    %% 抽取前差分时延
K;          %% 级联级数
fs;         %% 采样频率

A = [1,zeros(size(1:N-1)),-1]/N;  %% CIC滤波器的分子
B = [1,-1];                       %% CIC滤波器的分母

delt_f = Fpass/200;
NFFT = ceil(fs/2/delt_f);
if(Fpass+5*delt_f>Fstop)
    delt_f = delt_f/2;
    NFFT = ceil(fs/2/delt_f);
end
[G,f] = freqz(A,B,NFFT,fs);       %% 单级CIC滤波器的幅频特性
GK = G.^K;                        %% K级CIC滤波器的幅频特性
GK(1) = 1;                        %% 滤波器幅频第一个值为1

figure;
plot(f,20*log10(abs(GK)));title('cic频响')

xlabel('f,Hz'),ylabel('Gain, dB')
grid on
%% 补偿滤波
GK_comp = 1./abs(GK);              %% 补偿滤波器的幅频，幅频响应在通带内是CIC滤波器的倒数

FpassNum = ceil(NFFT*Fpass/fs*2);  %% 补偿滤波器通带频率对应的位置
if(mod(FpassNum,2) == 1)
    FpassNum = FpassNum+1;
end
fpass1 = f(2:FpassNum)';           %% 补偿滤波器通带频率
fstop1 = [Fstop fs/D/2];         %% 补偿滤波器阻带截止频率
dev = [(10^(Aps/20)-1)/(10^(Aps/20)+1) 10^(-Asp/20)];
Nord = firpmord([f(FpassNum),Fstop],[GK_comp(FpassNum),0],dev,fs/D);
% Nord = 54;                %% 补偿滤波器阶数
fo = [0 fpass1 fstop1];            %% 补偿滤波器总的频率

aopass = GK_comp(2:FpassNum)';    %% 补偿滤波器通带频率响应
ao = [1 aopass 0 0];               %% 补偿滤波器总的频率响应

cic_comp_flt1 = firpm(Nord,fo/(fs/D/2),ao);
cic_comp_flt = upsample(cic_comp_flt1,D);
[G_comp,f] = freqz(cic_comp_flt,1,NFFT,fs);

hold on
plot(f,20*log10(abs(G_comp)),'g');
hold on
plot(f,20*log10(abs(G_comp.*GK)),'r');
legend('CIC','Comp','CIC*Comp');
txt=['Delay=' num2str(M)  '; Decimate=' num2str(D)  '; Casc=' num2str(K)];

[G_comp1,f1] = freqz(cic_comp_flt1,1,NFFT,fs/D);

figure;
plot(f1,20*log10(abs(G_comp1)))
grid on


end