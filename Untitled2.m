%%
% 第一级
N = length(signal);% 数据点数
cic_sum1 = zeros(1,N);
cic_sum1(1) = signal(1);
for i = 2:N
    cic_sum1(i) = cic_sum1(i-1)+signal(i);
end
% 下抽取
cic_down = cic_sum1(1:5:end);
% 差分
cic_diff1 = zeros(1,N/5);
cic_diff1(1) = 0;
for i = 2:N/5
    cic_diff1(i) = cic_down(i) - cic_down(i-D);
end

% 第二级
cic_sum2 = zeros(1,N/25);
cic_sum2(1) = cic_diff1(1);
for i = 2:N/25
    cic_sum2(i) = cic_sum2(i-1)+cic_diff1(i);
end
% 下抽取
cic_down = cic_sum2(1:5:end);
% 差分
cic_diff2 = zeros(1,N/25);
cic_diff2(1) = 0;
for i = 2:N/R/R*5*5
    cic_diff2(i) = cic_down(i) - cic_down(i-D);
end

% 第三级
cic_sum3 = zeros(1,N/R/R*5*5);
cic_sum3(1) = cic_diff2(1);
for i = 2:N/R/R*5*5
    cic_sum3(i) = cic_sum3(i-1)+cic_diff2(i);
end
% 下抽取
cic_down = cic_sum3(1:R/5:end);
% 差分
cic_diff3 = zeros(1,N/R/R/R*5*5*5);
cic_diff3(1) = 0;
for i = 2:N/R/R/R*5*5*5
    cic_diff3(i) = cic_down(i) - cic_down(i-D);
end

% 第四级
cic_sum4 = zeros(1,N/R/R/R*5*5*5);
cic_sum4(1) = cic_diff3(1);
for i = 2:N/R/R/R*5*5*5
    cic_sum4(i) = cic_sum4(i-1)+cic_diff3(i);
end
% 下抽取
cic_down = cic_sum4(1:R/5:end);
% 差分
cic_diff4 = zeros(1,N/R/R/R/R*5*5*5*5);
cic_diff4(1) = 0;
for i = 2:N/R/R/R/R*5*5*5*5
    cic_diff4(i) = cic_down(i) - cic_down(i-D);
end

% 第五级
cic_sum5 = zeros(1,N/R/R/R/R*5*5*5*5);
cic_sum5(1) = cic_diff4(1);
for i = 2:N/R/R/R/R*5*5*5*5
    cic_sum5(i) = cic_sum5(i-1)+cic_diff4(i);
end
% 下抽取
cic_down = cic_sum5(1:R/5:end);
% 差分
cic_diff5 = zeros(1,N/R/R/R/R/R*5*5*5*5*5);
cic_diff5(1) = 0;
for i = 2:N/R/R/R/R/R*5*5*5*5*5
    cic_diff5(i) = cic_down(i) - cic_down(i-D);
end

out = cic_diff5;