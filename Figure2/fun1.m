function dx=fun1(t,x)

global b N Ki

dx=zeros(N,1);

%%% pulse pertubation
bB=0.5; bG=2;
if t>20 && t<60
    B=b;
    B(1)=bB;    B(3)=bG;
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
% ====