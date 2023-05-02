clc;clear all;close all;

G=zpk([],[-1 -3], [10])
Tm=(0.12)

Gd=c2d(G,Tm,'zoh')
Gd1=c2d(G,10*Tm,'zoh')

kp1=dcgain(G)
F1=feedback(G,1)
kp=dcgain(Gd)
F=feedback(Gd,1)

t=0:Tm:100*Tm %genera rampa

figure(6);grid on;hold on;lsim(F,t,t),title('Verificacion de error con entrada rampa');hold on

figure(7);
step(G);hold on;step(Gd);hold on;step(Gd1);hold on;grid on;

figure(1);
subplot(1,3,1),pzmap(G),title('Ceros y Polos - Continua G')
subplot(1,3,2),pzmap(Gd),title('Ceros y Polos - Discreta Gd')
subplot(1,3,3),pzmap(Gd1),title('Ceros y Polos - Discreta Gd1')

figure(2);
subplot(1,3,1),step(G),grid on,title('Rta escalon - Continua G')
subplot(1,3,2),step(Gd),grid on,title('Rta escalon - Discreta Gd')
subplot(1,3,3),step(Gd1),grid on,title('Rta escalon - Discreta Gd1')

figure(3);
step(F1);hold on;step(F);title('Rta escalon lazo cerrado')

figure(4);
subplot(1,2,1),rlocus(G),title('Rta escalon - Continua G')
subplot(1,2,2),rlocus(Gd),title('Rta escalón - Discreta Gd')

figure(5)
subplot(1,2,1),rlocus(Gd),title('Rta escalon - Discreta Gd')
subplot(1,2,2),rlocus(Gd1),title('Rta escalon - Discreta Tm*10')

%%
clc;clear all;close all;

G=zpk([],[-1 -1], [10])
Tm=(0.07)

% Gd=c2d(G,Tm,'zoh')
Gd=c2d(G,10*Tm,'zoh')
% 
kp1=dcgain(G)
F1=feedback(G,1)
kp=dcgain(Gd)
F=feedback(Gd,1)

t=0:Tm:100*Tm %genera rampa

% figure(6);grid on;hold on;lsim(F,t,t),title('Verificacion de error con entrada rampa');hold on
% 
figure(7);
step(F);hold on;grid on;

figure(8);
step(F1);hold on;grid on;
% 
% figure(1);
% subplot(1,3,1),pzmap(G),title('Ceros y Polos - Continua G')
% subplot(1,3,2),pzmap(Gd),title('Ceros y Polos - Discreta Gd')
% subplot(1,3,3),pzmap(Gd1),title('Ceros y Polos - Discreta Gd1')
% % 
% figure(2);
% subplot(1,3,1),step(G),grid on,title('Rta escalon - Continua G')
% subplot(1,3,2),step(Gd),grid on,title('Rta escalon - Discreta Gd')
% subplot(1,3,3),step(Gd1),grid on,title('Rta escalon - Discreta Gd1')

% figure(3);
% step(F1);hold on;step(F);title('Rta escalon lazo cerrado')
% 
% figure(4);
% subplot(1,2,1),rlocus(G),title('Rta escalon - Continua G')
% subplot(1,2,2),rlocus(Gd),title('Rta escalón - Discreta Gd')
% % 
% figure(5);
% subplot(1,2,1),rlocus(Gd),title('Rta escalon - Discreta Gd')
% subplot(1,2,2),rlocus(Gd1),title('Rta escalon - Discreta Tm*10')
