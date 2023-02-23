function [ odata ] = jiechanpi(idata)

num = length(idata);
ddata = diff(idata);
odata(1,1) = idata(1);
n = 0;

for i = 2:num  
    
    if(ddata(i-1)>pi)       
        n = n-1;
    elseif(ddata(i-1)<-pi)
        n = n+1;
    end
    
    odata(i,1)=idata(i)+n*2*pi;
    
end

end

