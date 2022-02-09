function dx=fun3(t,x)
global bb  N Ki

dx=zeros(N,1);
%%% OUP pertubation
    if t==0        
    B=bb(1,1:N);        
    else
    B=bb(ceil(t),1:N);        
    end
%%%

for i=1:N
dx(i)=x(i)*(B(i).*fi_Xk(i, x)-Ki(i).*x(i));
end
end


function fi=fi_Xk(i, x)
global n N Kij
fi=1;
K=1:N;K(i)=[];
for j=1:N-1
    k=K(j);
fi=fi*(Kij(i,k).^n/(Kij(i,k).^n+x(k).^n));
end
end
