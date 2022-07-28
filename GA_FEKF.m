%Guel-Cortez et al 2022. 
%Parameter Estimation of Fractional-Order Systems via Evolutionary Algorithms and the Extended Fractional Kalman Filter

close all;
clearvars;
clc

timerVal=tic;

L = 10;
Tspan=50;
T=0.05;
t=0:T:Tspan;
ks=2;
b=0.5;
N=[.8;.8];
Q=1e-5*eye(2,2);
R=1e-3;
u=sin(t);
[t,y]=FOsystem(ks,b,N,Q,R,T,t,L,u);
Q=1e-6*eye(5,5);
realp=[ks,b,N(1)];

options = optimoptions('ga','PlotFcn',"gaplotbestf",'UseParallel', true,'MaxStallGenerations',500,'MaxGenerations',1000);
[x,fval] = ga(@(gains) FOKFilter(t,u,y,gains,Q,R,T,L),3,[],[],[],[],[0;0;0],[3;1;1],[],options);

[RMSE,te,ye,x_e]=FOKFilter(t,u,y,x,Q,R,T,L);
figure 
set(gcf,'color','w');
plot(t,y,te,ye)
xlabel('t (s)','Interpreter','Latex','FontSize', 12)
ylabel('$y(t)/x_e(t)$','Interpreter','Latex','FontSize', 12)
legend('Measurement','Response with estimated parameters')
% 
 figure
   set(gcf,'color','w');
subplot(3,2,[1 2])
plot(t,y(1,:),t,x_e(1,:))
xlabel('t (s)','Interpreter','Latex','FontSize', 12)
ylabel('$y(t)/x_e(t)$','Interpreter','Latex','FontSize', 12)
legend('Measurement','Estimated')
subplot(3,2,3)
plot(t,x_e(2,:))
xlabel('t (s)','Interpreter','Latex','FontSize', 12)
ylabel('$\dot{x}_e(t)$','Interpreter','Latex','FontSize', 12)
subplot(3,2,4)
plot(t,x_e(3,:))
xlabel('t (s)','Interpreter','Latex','FontSize', 12)
ylabel('$k_e(t)$','Interpreter','Latex','FontSize', 12)
subplot(3,2,5)
plot(t,x_e(4,:))
xlabel('t (s)','Interpreter','Latex','FontSize', 12)
ylabel('$b_e(t)$','Interpreter','Latex','FontSize', 12)
subplot(3,2,6)
plot(t,x_e(5,:))
%ylim([N(1)-0.1,N(1)+0.1])
xlabel('t (s)','Interpreter','Latex','FontSize', 12)
ylabel('$\alpha_e(t)$','Interpreter','Latex','FontSize', 12)

estimatedp=[x_e(3:5,end)];
norm(realp-estimatedp,2)
Elapsetime=toc(timerVal);
