%% Figure 5 panel c (Ternary)
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

clear
clc
global n N Ki

%% Input 
% Coefficients and Conditions

%number of samples
Number=50; 

% Hill coefficient
n=2; 

% number of species
N=15; 

% dividing to 3 groups if N is a multiplie of 3
N3=N/3; 
Blue=1:N3;Red=N3+1:2*N3;Green=2*N3+1:3*N3;

% cells of interaction matrix
Kin=1; % intra-group inhibition
Kout=.6; % inter-group inhibition

% order of derivatives
mu=ones(N,1);
mu(Blue)=.6; % 1-Memory_B; For other species: mu(Red)=1-Memory_R, mu(Green)=1-Memory_G;

% death rate
Ki=1*ones(N,1); 

% t0=start time, T=the final time
t0=0; T=300; 

% time step size for computing
h=0.01; 

%% Computation

[b,r,g]=method3(mu,t0,T,h,Number,Kout,Kin,Blue,Red,Green);

%% plotting

figure      

% plot of the triangle
xx=[0  1 1/2 0];
yy=[0  0 sqrt(3)/2 0];
p1=plot(xx,yy,'k');
set(p1,'LineWidth',1.3)

% Find the location of the steady state for each sample and plot it as a colored dot
hold on
for i=1:Number
    G=g(i);B=b(i);R=r(i);
if G>B && G>R
pg=plot(1/2*(2*R+G),sqrt(3)/2*(G),'g.');
set(pg,'MarkerSize',10)
elseif R>B && R>G
pr=plot(1/2*(2*R+G),sqrt(3)/2*(G),'r.');
set(pr,'MarkerSize',10)
elseif B>G && B>R
pb=plot(1/2*(2*R+G),sqrt(3)/2*(G),'b.');
set(pb,'MarkerSize',10)
end
end

axis([-.2,1.2,-.2,1.2])

set(gca,'XColor','none','YColor','none') % removing axes lines

% Text which species on plot
text(-.07,-.03,'X_{B}','FontSize',16)
text(1.024,-.02,'X_{R}','FontSize',16)
text(.4,.95,'X_{G}','FontSize',16)
hold off

title(['Memory_B=',num2str(1-mu(1)),...
    ',  Memory_R=',num2str(1-mu(6)),...
    ',  Memory_G=',num2str(1-mu(15))]);
