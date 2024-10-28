function [ result ] = cqj_BsplineBasisFunction(i,k,u,t)
if(k==0)
    if(t(i+1)<=u && u<t(i+2)) || (u==1 && t(i+2)==1)%注意1=t(i)<=u<t(i+1)=1时的情况,这里要用t<=u(i+1);
        result=1;
        return;
    else
        result=0;
        return;
    end
else
    if(t(i+k+1)-t(i+1)==0)
        alpha=0;
    else
        alpha=(u-t(i+1))/(t(i+k+1)-t(i+1));
    end
    if(t(i+k+2)-t(i+2)==0)
         beta=0;
    else
        beta=(t(i+k+2)-u)/(t(i+k+2)-t(i+2));
    end
end
 result=alpha*cqj_BsplineBasisFunction(i,k-1,u,t)+beta*cqj_BsplineBasisFunction(i+1,k-1,u,t);