function D = braycd(X1,X2)

% input : X1 and X2 as two vectors containing the abundances of species in Time 1 and Time 2
% output : the bray curtis distance
%         braycd= sigma(|X1-X2|)/(sigma(X1)+sigma(X2))


X=abs(X1-X2);
Cij=sum(X);

Si_plus_Sj=sum(X1)+sum(X2);

D=Cij/Si_plus_Sj;