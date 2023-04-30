
clc;close all;
G=zpk([],[-1 -3],[10])
t=0.12
Gd=c2d(G,t,'zoh')
T=0.12
% sisotool(Gd)
C=zpk([0.564 0.807],[1 0],0.93928,t) %muestra el compensador importado de sisotool
F=feedback(C*Gd,1) % sistema de lazo cerrado
pole(F)
zero(F)
pzmap(F)

figure(1)
step(F);hold on;
% plot(tout,yout(:,1))
figure(3)
plot(tout,yout(:,4:6));hold on;legend('Derivariva','Integral','Prop');