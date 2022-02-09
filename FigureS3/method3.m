function [BB,RR,GG]=method3(mu,t0,T,h,Number,Kout,Kin,Blue,Red,Green)
% -------------------------------------------------------------------------
%                   This code several times solves an N species microbial community model
%                   described by fractional differential equations:
%                   D^mu(Xi)=Xi(bi.Fi-ki.Xi)
%                   where Fi=\prod[Kik^n/(Kik^n+Xk^n)], k=1,...,N and k~=i
%                   D is the fractional Caputo derivative and mu is its order                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs               
%        Number - Number of times the problem will be solved; e.g. Number=50
%        ------------------------------------------------------------------
%        mu - Order of derivatives, e.g. mu=ones(1,N);  % 1-Memory
%                                        Or
%                                        mu(Blue)=0.6; %  1-Memory_B
%                                        mu(Red)=1;      %  1-Memory_R
%                                        mu(Green)=1;    %  1-Memory_G 
%        ------------------------------------------------------------------
%        n -  Hill coefficient, e.g. n=2;                              
%        ------------------------------------------------------------------
%        N -  Number of Species (N=3k for specifying to 3 groups), e.g. N=15;
%        ------------------------------------------------------------------    
%        Kin - intra-group inhibition for Matrix interaction Kij e.g. Kin=1;
%        Kout - inter-group inhibition for Matrix interaction Kij, e.g. Kout=0.6;
%        ------------------------------------------------------------------
%        Ki - Death rate, e.g. Ki=1*ones(N,1); 
%        ------------------------------------------------------------------
%        T - Final time, e.g. T=700;
%        t0 - Initial time, e.g. t0=0;
%        ------------------------------------------------------------------
%        h - time step size for computing, e.g. h=0.01; 
%--------------------------------------------------------------------------
%
%*NOTE* No need to define the grwoth rates (b) and initial values (x0). 
%       They are generated randomly for each sample in 'method3' by:
%       b=ones(N,1)+.1*randn(N,1); % growth rate
%       x0=.1*rand(N,1); % initial conditions
%-----------------------------------
% Outputs
%        Ternary plot showing steady state distribution of several samples
%-----------------------------------
% Moein Khalighi - September 2020
%
%
%  Please, report any problem or comment to :
%          moein dot khalighi at utu dot fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global N b Kij

BB=zeros(Number,1);RR=zeros(Number,1);GG=zeros(Number,1);

for num=1:Number
tStart=tic;
    
% interaction matrix
var=.1*randn(N);
Kij=Kout*ones(N);
Kij(Blue,Blue)=Kin; Kij(Red,Red)=Kin; Kij(Green,Green)=Kin; % same coexistance
Kij=Kij+var;

b=ones(N,1)+.1*randn(N,1); % growth rate
x0=.1*rand(N,1); % initial conditions

% solver for fractional differential equation
[~, x] = FDE_PI12_PC(mu,@fun,t0,T,x0,h); 

% Extract values for the location of steady states in the triangle
G=sum(x(Green,end)); 
B=sum(x(Blue,end)); 
R=sum(x(Red,end)); 
all=G+R+B;
GG(num)=G/all;RR(num)=R/all;BB(num)=B/all;

tEnd=toc(tStart);
disp(['Computation time for sample (', num2str(num), ') is ',num2str(tEnd)])
end
end
% =========================================================================
% =========================================================================
function dx=fun(~,x)

global b N Ki

dx=zeros(N,1);

for i=1:N
dx(i)=x(i)*(b(i).*fi_Xk(i, x)-Ki(i).*x(i));
end
end
% =========================================================================
% =========================================================================
function fi=fi_Xk(i, x)
global n N Kij
fi=1;
K=1:N;K(i)=[];
for j=1:N-1
    k=K(j);
fi=fi*(Kij(i,k).^n/(Kij(i,k).^n+x(k).^n));
end
end