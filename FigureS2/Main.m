%% Figure S2
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

%% Inputs
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

% death rate
Ki=1*ones(N,1); 

% t0=start time, T=the final time
t0=0; T=300; 

% time step size for computing
h=0.01; 

%% Computation

% order of derivatives
mu=ones(N,1);

ind=[1;0.9;0.6];
Indx=permn(ind,3);
figure

for J=1:length(Indx)
mu(Blue)=Indx(J,1); % 1-Memory_B; For other species: mu(Red)=1-Memory_R, mu(Green)=1-Memory_G;
mu(Red)=Indx(J,2); % 1-Memory_B; For other species: mu(Red)=1-Memory_R, mu(Green)=1-Memory_G;
mu(Green)=Indx(J,3); % 1-Memory_B; For other species: mu(Red)=1-Memory_R, mu(Green)=1-Memory_G;

[b,r,g]=method3(mu,t0,T,h,Number,Kout,Kin,Blue,Red,Green);


%%plotting

%% pie chart

subplot(9,9,3*J-1:3*J)

PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];
ax = gca(); 
pieData = [.3 .4 .3]; 
hh=pie(ax,[mean(b), mean(r), mean(g)]);

newColors = [[0,0,1];PcR;PcG];

patchHand = findobj(hh, 'Type', 'Patch'); 
% Set the color of all patches using the nx3 newColors matrix
set(patchHand, {'FaceColor'}, mat2cell(newColors, ones(size(newColors,1),1), 3))
% Or set the color of a single wedge
patchHand(1).FaceColor = 'b';
patchHand(2).FaceColor = PcR;
patchHand(3).FaceColor = PcG;

% title('Propotion of abundance', 'fontsize',18)

set(findobj(hh,'type','text'),'fontsize',12);
%% Bar plots


subplot(9,9,3*J-2)

for i = 1:3    
            if i==1
                            p1(i)=bar(i,1-mu(1));
                            hold on
                set(p1(i),'FaceColor','b');
            elseif i==2         
                            p1(i)=bar(i,1-mu(6));
                set(p1(i),'FaceColor',PcR);
            elseif i==3
                            p1(i)=bar(i,1-mu(15));
                set(p1(i),'FaceColor',PcG);                
            end
end
  axis([.5 3+0.5 0 round(max(mu),1)])% adjusting axis
      hold off      
%                       ylabel('{Memory}','FontSize',18)
                      set(gca,'FontSize',18)
                      set(gca,'xticklabel',{})
                      
                      %%
end
