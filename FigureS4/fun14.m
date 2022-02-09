function dx=fun14(t,x)

global b N Ki

dx=zeros(N,1);

%%% pulse pertubation
b11=0.2; b12=4.5;
if t>60 && t<100
    B=b;
    B(1)=b11;
elseif t>400 && t<530
    B=b;
    B(1)=b12;
else
    B=b;
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