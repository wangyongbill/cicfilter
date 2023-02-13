function Am = beisaier(x,x0,Jn,ifplot)

% 输出Jn阶贝塞尔函数图像,包含0阶的贝塞尔
% 以及特定位置点的各阶贝塞尔值
% 输入：
%     x  为函数图像长度
%     x0 为所需点的横坐标
% 输出：
%      Am 为特定点的各阶贝塞尔值

  len_x = length(x);
  Am_Jn = zeros(Jn+1,len_x);
  Am_x0 = zeros(Jn+1,1);
  
  for i = 1:Jn+1
      Am_Jn(i,:) = besselj(i-1,x);  % Jn阶贝塞尔函数值     
  end 
    
  for i = 1:Jn+1
      Am_x0(i,:) = besselj(i-1,x0); % Jn阶贝塞尔函数值     
  end
  Am = Am_x0;
  
  if 1== ifplot
      figure;plot(x,Am_Jn);grid on;title('bessel函数');hold on;
      plot([x0 x0],[-0.5 1],'Marker','none');hold on
      plot(x0,Am_x0,'o');hold on
      text(x0,-0.5,num2str(x0));
      for i = 1:Jn+1
          text(x0+0.01,Am_x0(i),['第' num2str(i-1) '阶: ' num2str(Am_x0(i))]);
      end
  end
 
end