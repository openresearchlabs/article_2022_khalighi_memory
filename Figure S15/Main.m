%% Figure S15
clear
clc

global n N b Ki Kij

%% inputs
X00=[2 .01 .01
    .01 2 .01
    .01 .01 2]; % initial conditions

n=2; % Hill coefficient
N=3; % number of species
Kij=0.1*ones(N); % interaction matrix
Ki=1*ones(N,1); % death rate
t0=0; T=1000; % t0 is start time, and T is the end time

Bb=0.01:.05:5;
% Bb=0.01:.1:1;

% alpha=[.6,.6,1]; % order of derivatives
alpha=[1,1,1]; % order of derivatives

h=0.01; % step size for computing
for j=1:length(Bb)
    tic
    b=[Bb(j), .95, 1.05]; % growth rate

for i=1:3
    x0=X00(i,:); % initial conditions
%     if Bb(j)>4.7 && i==2
%             x0(i)=1;
% [t, x] = FDE_PI12_PC(alpha,@fun,t0,800,x0',h,[],2);        
%     elseif Bb(j)>4.7 && i==3
%             x0(i)=1;
% [t, x] = FDE_PI12_PC(alpha,@fun,t0,800,x0',h,[],2);        
    if Bb(j)<1 && i==1
            x0(i)=.2;
[t, x] = FDE_PI12_PC(alpha,@fun,t0,1000,x0',h,[],2);
    else
[t, x] = FDE_PI12_PC(alpha,@fun,t0,T,x0',h,[],2);
    end
RelX=x./(ones(N,1)*sum(x)); % Relative abundances

XB(i,j,:)=RelX(:,end);
end
    toc
end

%% plotting

f=figure;
f.Position = [0 0 1800 500];
f.Renderer='painters';
% defining blindfriendly colors (red and green)
PcR= [0.92,0.27,0.18];
PcG= [0.18,0.40,0.14];
   
%        ------------------------------------------------------------------    
    % Plot the  abundance of species
    hold on
    Xb(:,:)=XB(:,:,1);
    Xr(:,:)=XB(:,:,2);
    Xg(:,:)=XB(:,:,3);
subplot(3,1,1)    
for j=1:i
%     p=semilogy(tt,Xb(j,:),'b');p.Color(4)=.3;
    p=plot(Bb, Xb(j,:)','.b');p.Color(4)=.3;
    
    hold on
    set(p,'LineWidth',3)
end

    hold on
%     p3=semilogy(tt(1:100:end),Xb(7,1:100:end),'ob','MarkerSize', 10);
%     p4=semilogy(tt(1:500:end),Xb(6,1:500:end),'sqb','MarkerSize', 10);
    ylabel('Abundance', 'FontSize',28)
%     
%     yticks([0.01 0.1 1])
% yticklabels({'-2','-1','0'})
% Settings for the general plot

xlabel('Growth rate of the blue species')
set(gca,'FontSize',20)
% set(gca, 'YScale', 'log')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,2)

hold on
    for j=1:i
                p1=plot(Bb,Xr(j,:)','.','color',PcR);p1.Color(4)=.3;
            %     p1=semilogy(tt,Xr(j,:),'color',PcR);p1.Color(4)=.3;
        
    set(p1,'LineWidth',3)
    end
%     semilogy(tt(1:100:end),Xr(7,1:100:end),'o','color',PcR,'MarkerSize', 10);
%     semilogy(tt(1:100:end),Xr(6,1:100:end),'sq','MarkerSize', 10,'LineWidth',2);
    
%     
%     yticks([0.01 0.1 1])
% yticklabels({'-2','-1','0'})
xlabel('Growth rate of the blue species')
set(gca,'FontSize',20)
box on
% set(gca, 'YScale', 'log')
%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(3,1,3)

hold on
    for j=1:i
       
%     p2=semilogy(tt,Xg(j,:),'color',PcG);p2.Color(4)=.3;
p2=plot(Bb,Xg(j,:)','.','color',PcG);p2.Color(4)=.3;
        
    set(p2,'LineWidth',3)
    end
%     semilogy(tt(1:100:end),Xg(7,1:100:end),'om','MarkerSize', 10);
%     semilogy(tt(1:100:end),Xg(6,1:100:end),'sqk','MarkerSize', 10,'LineWidth',2);
    
    
%     yticks([0.01 0.1 1])
% yticklabels({'-2','-1','0'})
xlabel('Growth rate of the blue species')
% axis([0,T,.001,1])
set(gca,'FontSize',20)
box on
% set(gca, 'YScale', 'log')