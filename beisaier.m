function Am = beisaier(x,x0,Jn,ifplot)

% ���Jn�ױ���������ͼ��,����0�׵ı�����
% �Լ��ض�λ�õ�ĸ��ױ�����ֵ
% ���룺
%     x  Ϊ����ͼ�񳤶�
%     x0 Ϊ�����ĺ�����
% �����
%      Am Ϊ�ض���ĸ��ױ�����ֵ

  len_x = length(x);
  Am_Jn = zeros(Jn+1,len_x);
  Am_x0 = zeros(Jn+1,1);
  
  for i = 1:Jn+1
      Am_Jn(i,:) = besselj(i-1,x);  % Jn�ױ���������ֵ     
  end 
    
  for i = 1:Jn+1
      Am_x0(i,:) = besselj(i-1,x0); % Jn�ױ���������ֵ     
  end
  Am = Am_x0;
  
  if 1== ifplot
      figure;plot(x,Am_Jn);grid on;title('bessel����');hold on;
      plot([x0 x0],[-0.5 1],'Marker','none');hold on
      plot(x0,Am_x0,'o');hold on
      text(x0,-0.5,num2str(x0));
      for i = 1:Jn+1
          text(x0+0.01,Am_x0(i),['��' num2str(i-1) '��: ' num2str(Am_x0(i))]);
      end
  end
 
end