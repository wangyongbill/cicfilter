%  Design a minimum-order CIC compensator that compensates...
%  for the droop in the passband for the CIC decimator.
Fs = 160e6;    % Input sampling frequency
Fpass = 3.2e6; % Frequency band of interest
D = 25;        % Decimation factor of CIC
d1 = fdesign.decimator(D,'CIC',1,Fpass,65,Fs); %design a cic filter
Hcic = design(d1);

Hd(1) = cascade(dfilt.scalar(1/gain(Hcic)),Hcic);
d2 = fdesign.ciccomp(Hcic.DifferentialDelay, ...
    Hcic.NumberOfSections,Fpass,1.625e6,.005,66,Fs/D); % design a cic compensator filter
Hd(2) = design(d2);
fcfwrite([Hcic Hd(2)],'CICdesciption','dec'); % 其中，生成的.fcf文件描述滤波器的结构
hvt=fvtool(Hd(1),Hd(2),cascade(Hd(1),Hd(2)),'Fs',[Fs Fs/D Fs], ...   % plot whole response
    'ShowReference', 'off');
legend(hvt, 'CIC','CIC compensator', 'Whole response','Location', 'Northeast');
