%% Figure 2 
%        ------------------------------------------------------------------
%                   This code solves a there species microbial community model
%                   described by fractional differential equations:
%                   D^mu(Xi)=X_i(bi.Fi-ki.Xi)
%                   where Fi=\prod[Kik^n/(Kik^n+Xk^n)], k=1,...,N and k~=i
%                   D is the fractional Caputo derivative and mu is its order  
%
%  For a article titled "Quantifying the impact of ecological memory on the dynamics of interacting communities"         
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
%        b - Growth rates (including pulse perturbations), e.g. b=[1, .95, 1.05];
%
%---------------------------------------
% Outputs
%        t - Simulated time interval
%        x - Species abundances 
%        B - Growth rates including perturbation
%---------------------------------------
%
%
%  Please, report any problem or comment to :
%          moein dot khalighi at utu dot fi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
clc
global n N Ki b Kij T x0

%% Inputs
% Coefficients and Conditions

% [mu_B,mu_R,mu_G] Order of derivatives,  0<mu(i)=<1
mu1=[1,1,1];    % No memory
mu2=.9*[1,1,1]; % With memory

n=2; % Hill coefficient

N=3; % number of species

Kij=0.1*ones(N); % interaction matrix

Ki=1*ones(N,1); % death rate

T=350; %  final time

b=[1, .95, 1.05]; % growth rates (Pulses are applied in fun1.m and fun12.m)

x0=[.99;.01;.01]; % initial conditions

t0=0; % initial time
h=0.01; % step size for computing

Fun1=@fun1;  % ODE model including a pulse for panel a 
Fun2=@fun12; % ODE model including a pulse for panel b
  
% solver for fractional differential equation
[t1, x1] = FDE_PI12_PC(mu1,Fun1,t0,T,x0,h); % results for panel a, no memory
[t2, x2] = FDE_PI12_PC(mu2,Fun1,t0,T,x0,h); % results for panel a, with memory
[~, x12] = FDE_PI12_PC(mu1,Fun2,t0,T,x0,h); % results for panel b, no memory
[~, x22] = FDE_PI12_PC(mu2,Fun2,t0,T,x0,h); % results for panel b, with memory
%% plotting

% defining blindfriendly colors (red and green)
PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];

%%     plotting

f1=figure;
f1.Renderer='painters';

tt=t1;
    
subplot(2,3,1)
%plot of growth rates
    ttb=0:.1:T;Nt=length(ttb); % time simulating with 0.1 step size
         
%     % Simulate growth rates with pulse
    BB=zeros(1,Nt);RR=zeros(1,Nt);GG=zeros(1,Nt);
    for i=1:Nt
        if ttb(i)>20 && ttb(i)<60    
        BB(i)=.5;GG(i)=2;
        else
        BB(i)=1;        GG(i)=b(3);
        end
        RR(i)=b(2);
    end
% %     


ff1 = [1 2 3 4];

    v1 = [20 0.5; 60 0.5; 60 2; 20 2];
         h2=patch('Faces',ff1,'Vertices',v1,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
        % Remove highlighted bars from the legend
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(gca,'YScale','log')
    semilogy(ttb,BB,'Color',[0,0,1],'LineWidth',4);
    hold on
    semilogy(ttb,RR,'Color',PcR,'LineWidth',4);
    semilogy(ttb,GG,'Color',PcG,'LineWidth',4);
%     Settings for plot
axis tight
    
    ylabel('Log growth rates')
    xlabel('Time')
    legend('b_B\approx1','b_R=0.95', 'b_G\approx1.05')
    set(gca,'FontSize',23)
    
     yticks([10^-.3 10^0 10^.3])
yticklabels({'-0.30','0','0.3'})
        set(gca,'YScale','log')
    
    Pos = [795 359 708 457];
set(0, 'DefaultFigurePosition', Pos);
    
    %%plotting relative abundance of species
    
    subplot(2,3,2)
    
    ff = [1 2 3 4];

    v = [20 .009; 60 .009; 60 1+.1; 20 1+.1];

   h12=patch('Faces',ff,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on


        set(gca,'YScale','log')
    p12=semilogy(tt,x1(1,:),'Color',[0,0,1]);
    hold on
    set(p12,'LineWidth',4)
    p23=semilogy(tt,x1(2,:),'Color',PcR);
    set(p23,'LineWidth',4)
    p1=semilogy(tt,x1(3,:),'Color',PcG);
    set(p1,'LineWidth',4)
    
Pos = [795 359 708 457];
set(0, 'DefaultFigurePosition', Pos);
    
    ylabel('Log abundance')
    
    yticks([10^-2 10^-1 10^0])
yticklabels({'-2','-1','0'})

    % Remove highlighted bars from the legend
%     set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% hleg1=legend('show');
% set(hleg1,'Location','bestoutside','FontSize',11)

set(gca,'FontSize',23)
xlabel('Time')
axis([0,T,0,1.03])

subplot(2,3,3)

    set(gca,'YScale','log')
% 
    p12=semilogy(tt,x1(1,:),'Color',[0,0,1]);p12.Color(4)=.3;
%     p12=semilogy(tt,X1(1,:),'Color',[0,0,1]);
    hold on
    set(p12,'LineWidth',4)
    p23=semilogy(tt,x1(2,:),'Color',PcR);p23.Color(4)=.3;
%     p23=semilogy(tt,X1(2,:),'Color',PcR);
    set(p23,'LineWidth',4)
    p1=semilogy(tt,x1(3,:),'Color',PcG);p1.Color(4)=.3;
%     p1=semilogy(tt,X1(3,:),'Color',PcG);
    set(p1,'LineWidth',4)
    
Pos = [795 359 708 457];
set(0, 'DefaultFigurePosition', Pos);
    
%     ylabel('Log abundance')
    
    yticks([10^-2 10^-1 10^0])
yticklabels({'-2','-1','0'})

set(gca,'FontSize',23)
xlabel('Time')
axis([0,T,0,1.03])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    h32=patch('Faces',ff,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
     set(gca,'YScale','log')
      p11=semilogy(tt,x2(1,:),'Color',[0,0,1]);
    hold on
    set(p11,'LineWidth',4)
    p22=semilogy(tt,x2(2,:),'Color',PcR);
   set(p22,'LineWidth',4)
    p=semilogy(tt,x2(3,:),'Color',PcG);set(p,'LineWidth',4);
    
    axis([0,T,0,1.03])
    set(gca,'FontSize',23)
xlabel('Time')

Pos = [795 359 708 457];
set(0, 'DefaultFigurePosition', Pos);

    yticks([10^-2 10^-1 10^0])
yticklabels({'-2','-1','0'})


    
subplot(2,3,4)
%plot of growth rates
    ttb=0:.1:T;Nt=length(ttb); % time simulating with 0.1 step size
%         
%     % Simulate growth rates with pulse
    BB=zeros(1,Nt);RR=zeros(1,Nt);GG=zeros(1,Nt);
    for i=1:Nt
        if ttb(i)>20 && ttb(i)<60    
        BB(i)=.5;GG(i)=2.2;
        else
        BB(i)=1;        GG(i)=b(3);
        end
        RR(i)=b(2);
    end
% %     


ff1 = [1 2 3 4];

    v1 = [20 0.5; 60 0.5; 60 2.2; 20 2.2];
         h2=patch('Faces',ff1,'Vertices',v1,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
        % Remove highlighted bars from the legend
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        set(gca,'YScale','log')
    pb1=semilogy(ttb,BB,'Color',[0,0,1],'LineWidth',4);
    hold on
    pb2=semilogy(ttb,RR,'Color',PcR,'LineWidth',4);
    pb3=semilogy(ttb,GG,'Color',PcG,'LineWidth',4);
%     Settings for plot
axis tight
    
    ylabel('Log growth rates')
    xlabel('Time')
%     legend('b_B\approx1','b_R=0.95', 'b_G\approx1.05')
    set(gca,'FontSize',23)
    
     yticks([10^-.3 10^0 10^.34])
yticklabels({'-0.30','0','0.34'})
        set(gca,'YScale','log')
    
    Pos = [795 359 708 457];
set(0, 'DefaultFigurePosition', Pos);
    
%%plotting relative abundance of species
    
    subplot(2,3,5)
    
    ff = [1 2 3 4];

    v = [20 .009; 60 .009; 60 1+.1; 20 1+.1];

   h1=patch('Faces',ff,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on


        set(gca,'YScale','log')
    p12=semilogy(tt,x12(1,:),'Color',[0,0,1]);
    hold on
    set(p12,'LineWidth',4)
    p23=semilogy(tt,x12(2,:),'Color',PcR);
    set(p23,'LineWidth',4)
    p1=semilogy(tt,x12(3,:),'Color',PcG);
    set(p1,'LineWidth',4)
    
Pos = [795 359 708 457];
set(0, 'DefaultFigurePosition', Pos);
    
    ylabel('Log abundance')
    
    yticks([10^-2 10^-1 10^0])
yticklabels({'-2','-1','0'})

set(gca,'FontSize',23)
xlabel('Time')
axis([0,T,0,1.03])

subplot(2,3,6)
% defining blindfriendly colors (red and green)


        hold on

    set(gca,'YScale','log')
% 
    p12=semilogy(tt,x12(1,:),'Color',[0,0,1]);p12.Color(4)=.3;
%     p12=semilogy(tt,X1(1,:),'Color',[0,0,1]);
    hold on
    set(p12,'LineWidth',4)
    p23=semilogy(tt,x12(2,:),'Color',PcR);p23.Color(4)=.3;
%     p23=semilogy(tt,X1(2,:),'Color',PcR);
    set(p23,'LineWidth',4)
    p1=semilogy(tt,x12(3,:),'Color',PcG);p1.Color(4)=.3;
%     p1=semilogy(tt,X1(3,:),'Color',PcG);
    set(p1,'LineWidth',4)
    
Pos = [795 359 708 457];
set(0, 'DefaultFigurePosition', Pos);
    
%     ylabel('Log abundance')
    
    yticks([10^-2 10^-1 10^0])
yticklabels({'-2','-1','0'})

set(gca,'FontSize',23)
xlabel('Time')
axis([0,T,0,1.03])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    h3=patch('Faces',ff,'Vertices',v,'FaceColor','k','EdgeColor','non', 'FaceAlpha',.16);
    hold on
     set(gca,'YScale','log')

     p11=semilogy(tt,x22(1,:),'Color',[0,0,1]);
    hold on
    set(p11,'LineWidth',4)

    p22=semilogy(tt,x22(2,:),'Color',PcR);
   set(p22,'LineWidth',4)

   p=semilogy(tt,x22(3,:),'Color',PcG);set(p,'LineWidth',4);
    
    axis([0,T,0,1.03])
    set(gca,'FontSize',23)
xlabel('Time')

Pos = [795 359 708 457];
set(0, 'DefaultFigurePosition', Pos);

    yticks([10^-2 10^-1 10^0])
yticklabels({'-2','-1','0'})

