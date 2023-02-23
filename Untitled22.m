intOut = zeros(N,1);
delay_intOut = zeros(N,1);
yn = [];
yn_tmp = [];
for i = 1:length(data_pdm)
    tmp = data_pdm(i);
    for j = 1:N      %%integrate
        intOut(j) = intOut(j) + tmp;
        tmp = intOut(j);
    end
    
    if mod(i,R) == 1  %%decimator
        tmp = intOut(N);
        yn_tmp = [yn_tmp tmp];
        for j = 1:N  %%comb
            combOut = tmp - delay_intOut(j,1);
            delay_intOut(j,1) = tmp;
            tmp = combOut ;
        end
        yn = [yn combOut];
    end
end
data_cic = yn(:)/(R^N);