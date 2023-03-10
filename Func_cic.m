function out = Func_cic(signal,M,R,D)
%%
Len = length(signal);       % 输入数据长度
MIntgrat = zeros(1,M);       % M级积分寄存器
MDiff = zeros(1,M);
IntgratOut = zeros(1,Len);  % 积分输出
DiffOut = zeros(1,Len);   % 差分输出
%% M级积分器，积分后位数展宽，
%  根据公式计算积分后的位数Bout = Bin-1+M*log2(R*D)
flg = 0;
for i = 1:Len
    indata = signal(i);
    for j = 1:M
        tmp = MIntgrat(j);
        MIntgrat(j) = indata + MIntgrat(j);
        indata = tmp;
    end
    IntgratOut(i) = MIntgrat(M);
end
%% 梳状滤波器 实现R抽取与D差分
% 抽取
DiffIn = IntgratOut(1:R:end);
for i = 1:length(DiffIn)
    difftemp = DiffIn(i);
    for j = 1:M
        difftemp2 = difftemp - MDiff(j);  % D = R差分时延等于抽取倍率
        MDiff(j) = difftemp;
        difftemp = difftemp2;
    end
    flg = flg+1;
    DiffOut(flg) = difftemp2;
end

out = DiffOut(1:flg);
%%
% N = length(signal);% 数据点数
%
% % 积分器级联
% cic_sum = zeros(1,N);
%
% cic_sum(1) = signal(1);
% for i = 2:N
%     cic_sum(i) = cic_sum(i-1)+signal(i);
% end
%
% % 下抽取
% cic_down = cic_sum(1:R:end);
%
% % 差分级联
% cic_diff = zeros(1,N/R);
%
% cic_diff(1) = 0;
% for i = 2:N/R
%     cic_diff(i) = cic_down(i) - cic_down(i-D);
% end
%
% out = cic_diff;
%% %%
% % 积分器级联
% cic_sum1 = zeros(1,N);
% cic_sum2 = zeros(1,N);
% cic_sum3 = zeros(1,N);
% cic_sum4 = zeros(1,N);
% cic_sum5 = zeros(1,N);
%
% cic_sum1(1) = signal(1);
% for i = 2:N
%     cic_sum1(i) = cic_sum1(i-1)+signal(i);
% end
%
% cic_sum2(1) = cic_sum1(1);%signal(1);
% for i = 2:N
%     cic_sum2(i) = cic_sum2(i-1)+cic_sum1(i);
% end
%
% cic_sum3(1) = cic_sum2(1);
% for i = 2:N
%     cic_sum3(i) = cic_sum3(i-1)+cic_sum2(i);
% end
%
% cic_sum4(1) = cic_sum3(1);
% for i = 2:N
%     cic_sum4(i) = cic_sum4(i-1)+cic_sum3(i);
% end
%
% cic_sum5(1) = cic_sum4(1);
% for i = 2:N
%     cic_sum5(i) = cic_sum5(i-1)+cic_sum4(i);
% end
%
% % 下抽取
% cic_down = cic_sum5(1:R:end);
%
% % 差分级联
% cic_diff1 = zeros(1,N/R);
% cic_diff2 = zeros(1,N/R);
% cic_diff3 = zeros(1,N/R);
% cic_diff4 = zeros(1,N/R);
% cic_diff5 = zeros(1,N/R);
% D = 2;
% cic_diff1(1) = 0;
% for i = D+1:N/R
%     cic_diff1(i) = cic_down(i) - cic_down(i-D);
% end
%
% cic_diff2(1) = 0;
% for i = D+1:N/R
%     cic_diff2(i) = cic_diff1(i) - cic_diff1(i-D);
% end
%
% cic_diff3(1) = 0;
% for i = D+1:N/R
%     cic_diff3(i) = cic_diff2(i) - cic_diff2(i-D);
% end
%
% cic_diff4(1) = 0;
% for i = D+1:N/R
%     cic_diff4(i) = cic_diff3(i) - cic_diff3(i-D);
% end
%
% cic_diff5(1) = 0;
% for i = D+1:N/R
%     cic_diff5(i) = cic_diff4(i) - cic_diff4(i-D);
% end
%
% out = cic_diff5;

end