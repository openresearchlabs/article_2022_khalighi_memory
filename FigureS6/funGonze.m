function dx=funGonze(t,x)

global b N Ki


dx=zeros(N,1);

for i=1:N
dx(i)=x(i)*(b(i).*fi_Xk(i, x)-Ki(i).*x(i));
end

end

