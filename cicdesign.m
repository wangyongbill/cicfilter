clear
clc
close all
%%
fs = 160e6;    % 采样频率
fpass = 75e5;  % 通带
fstop = 150e5;
R = 4;         % 抽取率
D = 1;         % 差分时延

% hcic = design(fdesign.decimator(R,'cic',D,fpass,13,fs),'SystemObject',true);
Dcic = fdesign.decimator(R,'cic',D,'Fp,Ast',fpass,60,fs);
hcic = design(Dcic,'SystemObject',true);
fvtool(hcic);


hd = design(fdesign.ciccomp(hcic.DifferentialDelay,hcic.NumSections,...
fpass,fstop,0.001,100,fs/R),'SystemObject',true);
fvtool(hd)

fvtool(hcic,hd,cascade(hcic,hd),'ShowReference','off','Fs',[fs fs/R fs])
legend('CIC Decimator','CIC Compensator','Resulting Cascade Filter');


