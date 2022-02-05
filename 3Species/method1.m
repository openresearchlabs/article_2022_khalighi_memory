function [t,x,B]=method1(Perturbation)
%        ------------------------------------------------------------------
%                   This code solves a there species microbial community model
%                   described by fractional differential equations:
%                   D^mu(Xi)=X_i(bi.Fi-ki.Xi)
%                   where Fi=\prod[Kik^n/(Kik^n+Xk^n)], k=1,...,N and k~=i
%                   D is the fractional Caputo derivative and mu is its order  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs                   
%        ------------------------------------------------------------------
%        mu - Order of derivatives, [mu_B,mu_R,mu_G]  0<mu(i)=<1, e.g. mu=[1,.2,1];
%        ------------------------------------------------------------------
%        n -  Hill coefficient, e.g. n=2;
%        ------------------------------------------------------------------
%        N -  Number of Species, e.g. N=3;
%        ------------------------------------------------------------------
%        Kij - Interation matrix, e.g. Kij=0.1*ones(N);
%        ------------------------------------------------------------------
%        Ki - Death rate, e.g. Ki=1*ones(N,1);
%        ------------------------------------------------------------------
%        T - Final time, e.g. T=600;
%        ------------------------------------------------------------------
%        x0 - Initial conditions, e.g. x0=[1/3;1/3;1/3];
%        ------------------------------------------------------------------
%        Perturbation - This is changes of the species growth rates, e.g. Perturbation='OUP';
%                       Possible usages 'False', 'Pulse, 'Periodic', 'OUP', and 'OUP_new'                      
%                       False: No perturbation
%                       Pulse1: A pulse in (20,60) for Figure 2a
%                       Pulse2: Similar to Pulse1 with a greater strength for Figure 2b
%                       Pulse3: Two pulses in (60,100) and (200,330) for Figure 3a
%                       Pulse4: Two pulses in (60,100) and (400,530) for Figure S2a
%                       Periodic: Periodic perturbation with 20 span for Figure 3b
%                       OUP: Stochastic pertubation used in the paper; requirement: T=<700, For Figure 4b-c & S2b                          
%                       OUP_new: New generating stochastic perturbation
%        ------------------------------------------------------------------
%        b - Growth rates for cases: False, Pulse, and Periodic, e.g. b=[1, .95, 1.05];
%
%---------------------------------------
% Outputs
%        t - Simulated time interval
%        x - Species abundances 
%        B - Growth rates including perturbation
%---------------------------------------
% Moein Khalighi - September 2020
%
%
%  Please, report any problem or comment to :
%          moein dot khalighi at utu dot fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global T N mu x0 bb b

switch Perturbation
    case 'False'
        Fun=@fun;
    %----------------------------------------------------------------------
    case 'Pulse1'
        Fun=@fun1;
    %---------------------------------------------------------------------- 
        case 'Pulse2'
        Fun=@fun12;
    %---------------------------------------------------------------------- 
        case 'Pulse3'
        Fun=@fun13;
        % Check if the final time passes the second pulse
        if T<330
                    warning('MATLAB:FinalTimeNotProper',...
                        'The final time T=%d should be more than 330 to pass the second pulse',T);
        end
    %---------------------------------------------------------------------- 
        case 'Pulse4'
        Fun=@fun14;
        % Check if the final time passes the second pulse
        if T<530
                    warning('MATLAB:FinalTimeNotProper',...
                        'The final time T=%d should be more than 530 to pass the second pulse',T);
        end
    %---------------------------------------------------------------------- 
    case 'Periodic'    
        Fun=@fun2;
    %----------------------------------------------------------------------        
    case 'OUP'        
        Fun=@fun3;
        clear b        
        load('b3OUP.mat');
        bb=@(t,N)b(t,N);
        % Check if the final time is suitable with produced OUP
        if T>700
            error('MATLAB:FinalTimeNotCompatible', ...
                'The final time T=%d should be less than 700', T);
        end
    %----------------------------------------------------------------------
    case 'OUP_new'  
        b=OUP(T,N);
        bb=@(t,N)b(t,N);
        Fun=@fun3;
end
B=b;

t0=0; % initial time
h=0.01; % step size for computing


% solver for fractional differential equation
[t, x] = FDE_PI12_PC(mu,Fun,t0,T,x0,h);
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
function dx=fun1(t,x)

global b N Ki

dx=zeros(N,1);

%%% pulse pertubation
bB=0.5; bG=2;
if t>20 && t<60
    B=b;
    B(1)=bB;    B(3)=bG;
else
    B=b;
end
%%%

for i=1:N
dx(i)=x(i)*(B(i).*fi_Xk(i, x)-Ki(i).*x(i));
end
end% =========================================================================
% =========================================================================
function dx=fun12(t,x)

global b N Ki

dx=zeros(N,1);

%%% pulse pertubation
bB=0.5; bG=2.2;
if t>20 && t<60
    B=b;
    B(1)=bB;    B(3)=bG;
else
    B=b;
end
%%%

for i=1:N
dx(i)=x(i)*(B(i).*fi_Xk(i, x)-Ki(i).*x(i));
end
end
% =========================================================================
% =========================================================================
function dx=fun13(t,x)

global b N Ki

dx=zeros(N,1);

%%% pulse pertubation
b11=0.2; b12=4.5;
if t>60 && t<100
    B=b;
    B(1)=b11;
elseif t>200 && t<330
    B=b;
    B(1)=b12;
else
    B=b;
end
%%%

for i=1:N
dx(i)=x(i)*(B(i).*fi_Xk(i, x)-Ki(i).*x(i));
end
end
% =========================================================================
% =========================================================================
function dx=fun14(t,x)

global b N Ki

dx=zeros(N,1);

%%% pulse pertubation
b11=0.2; b12=4.5;
if t>60 && t<100
    B=b;
    B(1)=b11;
elseif t>400 && t<530
    B=b;
    B(1)=b12;
else
    B=b;
end
%%%

for i=1:N
dx(i)=x(i)*(B(i).*fi_Xk(i, x)-Ki(i).*x(i));
end
end
% =========================================================================
% =========================================================================
function dx=fun2(t,x)

global b N Ki

dx=zeros(N,1);

%%% periodic pertubation
mm=20;m=ceil(mod(t/(mm*4),mm));
b11=0.2; b12=4.5;
if t==0
    B=b;
elseif t>=mm*(4*m-3) && t<mm*(4*m-2)
    B=b;
    B(1)=b11;
    elseif t>=mm*(4*m-1) && t<mm*(4*m)
    B=b;
    B(1)=b12;
else 
        B=b;
end
%%%

for i=1:N
dx(i)=x(i)*(B(i).*fi_Xk(i, x)-Ki(i).*x(i));
end
end
% =========================================================================
% =========================================================================
function dx=fun3(t,x)
global bb  N Ki

dx=zeros(N,1);
%%% OUP pertubation
    if t==0        
    B=bb(1,1:N);        
    else
    B=bb(ceil(t),1:N);        
    end
%%%

for i=1:N
dx(i)=x(i)*(B(i).*fi_Xk(i, x)-Ki(i).*x(i));
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
% =========================================================================
% =========================================================================
function b=OUP(T,N)
% Ornstein Uhlenbeck stochastic process
% $dX_t = \theta (\mu - X_t)dt + \sigma dW_t$

%Seed for RNG
vector=1:2^20;
seed=vector(randi(length(vector)));

% OU parameters:
%mu, s0, vol, theta
phi = 1;         % Reverting level
b0 = 1;           % Initial value
Vol = .1;   % volatility of each step
theta = .08;      % Speed of mean reversion
             
% Using simple discretisation
b = OrnsteinUhlenbeck(T, N, seed, phi, b0, Vol, theta);
 end
