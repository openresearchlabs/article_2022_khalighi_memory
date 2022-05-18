function dx=Jfun(t,x,p)

global A mu

%%%pertubation
if t>50 && t<80
    B=mu;
B(1)=p+mu(1);
else
    B=mu;
end
%%%

dx11= (B(1)+A(1,1)*x(1)+A(1,2)*x(2))+x(1)*(A(1,1));
dx12= x(1)*(A(1,2));
dx21= x(2)*(A(2,1));
dx22= (B(2)+A(2,1)*x(1)+A(2,2)*x(2))+x(2)*(A(2,2));

dx=[dx11 dx12;dx21 dx22];
end