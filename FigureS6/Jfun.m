function dx=Jfun(t,x)

global A mu

dx11= (mu(1)+A(1,1)*x(1)+A(1,2)*x(2))+x(1)*(A(1,1));
dx12= x(1)*(A(1,2));
dx21= x(2)*(A(2,1));
dx22= (mu(2)+A(2,1)*x(1)+A(2,2)*x(2))+x(2)*(A(2,2));

dx=[dx11 dx12;dx21 dx22];
end