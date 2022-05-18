function fi=fi_Xk(i, x)
global n N Kij
fi=1;
K=1:N;K(i)=[];
for j=1:N-1
    k=K(j);
fi=fi*(Kij(i,k).^n/(Kij(i,k).^n+x(k).^n));
end
end