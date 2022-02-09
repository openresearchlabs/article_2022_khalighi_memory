function dx=funGonze(t,x,p)

global b N Ki


%%%pertubation
if t>50 && t<100
    B=b;
B(1)=p+b(1);
else
    B=b;
end
%%%


dx=zeros(N,1);

for i=1:N
dx(i)=x(i)*(B(i).*fi_Xk(i, x)-Ki(i).*x(i));
end

end

