function dx=fun(t,x)

global A mu

dx1= x(1)*(mu(1)+A(1,1)*x(1)+A(1,2)*x(2));
dx2= x(2)*(mu(2)+A(2,1)*x(1)+A(2,2)*x(2));

dx=[dx1;dx2];
end

