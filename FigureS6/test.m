%% BU BT
clear 
clc
%%
global A mu


order1=[1,1];

ExpRelBT=[0.5873 0.62714 0.63223 0.66755 0.67419 0.70986];
ExpRelBU=1-ExpRelBT;
ExpT=0:12:12*length(ExpRelBT)-12;

% initial conditions
AbsIniBT=.02*ExpRelBT(1);
AbsIniBU=.02*ExpRelBU(1);
X0= [AbsIniBU; AbsIniBT];

mu=[0.599 0.626];

t0=0;
T=2000;
h=.1;
F=@fun;
JF=@Jfun;

A=[-0.9059 -0.9377;-0.972 -0.9597];

% Xx1=[(A(1,2)*mu(2)-A(2,2)*mu(1))/(A(1,1)*A(2,2)-A(1,2)*A(2,1)),...
%     (A(2,1)*mu(1)-A(1,1)*mu(2))/(A(1,1)*A(2,2)-A(1,2)*A(2,1))];
Xx1=[0,-mu(2)/A(2,2)];
Xx2=flip(Xx1);

[t1,x1]=FDE_PI2_IM(order1,F,JF,t0,T,X0,h);
% [t2,x2]=FDE_PI2_IM(order2,F,JF,t0,T,X0,h);

ER1= braycd(x1,Xx1');
% ER2= braycd(x2,Xx2');
%%
figure
% 
X=x1./sum(x1);
plot(t1(1:20:end),X(1,1:20:end),'b',t1(1:20:end),X(2,1:20:end),'r','LineWidth',8)
% hold on
% plot(t2,x2(1,:),'--b',t2,x2(2,:),'--r')

xlabel('Time (hr)')
ylabel('Abundance')
% legend('X1; memory = 0','X2; memory = 0')
set(gca,'FontSize',30)

% plot([t1(1) t1(end)],[Xx1(1) Xx1(1)],'k')
% plot([t1(1) t1(end)],[Xx1(2) Xx1(2)],'k')
% figure
% semilogy(t2,ER1,'--b',t2,ER2,'--r')


