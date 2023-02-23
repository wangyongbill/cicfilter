function [Out_data,y] = CIC_decimate_v1(In_Data,N,K)

    % �������ֲ��� 
    y_integral = zeros(K,2);
    for i = 1:length(In_Data)
        y_integral(1,1) =  In_Data(i) + y_integral(1,2);
        for j = 2:K
            y_integral(j,1) =  y_integral(j-1,2) + y_integral(j,1); % �������ֵ�Ԫ
        end
        y_integral(:,2) = y_integral(:,1); % ��ʱ��Ԫ
        
        y(i) = y_integral(K,2);  % ���ֲ������
    end
    
    % ������״����
    y = y(1:N:end);   % ��ȡ
    y_comb = zeros(K,2);
    Out_data = zeros(1,length(y));
    for i = 2:length(y)
        y_comb(1,2) = y(i-1); 
        y_comb(1,1) = y(i); 
        for j = 1:K-1
            y_comb(j+1,1) = y_comb(j,1) - y_comb(j,2);
        end
        Out_data(i) = y_comb(K,1) - y_comb(K,2);  
        y_comb(:,2) = y_comb(:,1);  % ��ʱ   
    end

end