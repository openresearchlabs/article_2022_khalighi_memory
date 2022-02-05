function dx=fun3(t,x)
global b N Ki T1 T2

%%%pertubation
if t>T1 && t<T2
    B=b;
%     B(indx)=b12;
%     B(indx)=.5*b(indx);
B(1:7)=3*b(1:7);
else
    B=b;
end
%%%

dx=zeros(N,1);

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

