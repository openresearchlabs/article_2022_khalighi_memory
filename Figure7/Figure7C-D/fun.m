function dx=fun(t,x,p)

global A mu

%%%pertubation
if t>50 && t<80
    B=mu;
B(1)=p+mu(1);
else
    B=mu;
end
%%%

dx1= x(1)*(B(1)+A(1,1)*x(1)+A(1,2)*x(2));
dx2= x(2)*(B(2)+A(2,1)*x(1)+A(2,2)*x(2));

dx=[dx1;dx2];
end

