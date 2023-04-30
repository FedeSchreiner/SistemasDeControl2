clc;clear all;close all;
syms x1 x1p x2 x2p theta thetap thetapp m b l g u e delta
x1 = e;
x2 = thetap;
x1p = thetap;
% x2p = (1/(m*l^2))*(u-b*x2-g*m*l*sin(x1+delta)); % no lineal
x2p = (1/(m*l^2))*(u-b*x2-g*m*l*(x1+delta)); % lineal

A = [[subs(subs(diff(x1p,x1),x1,0),x2,0) subs(subs(diff(x1p,x2),x1,0),x2,0)];
     [subs(subs(diff(x2p,x1),x1,0),x2,0) subs(subs(diff(x2p,x2),x1,0),x2,0)]];
B = [subs(subs(diff(x1p,u),x1,0),x2,0);
     subs(subs(diff(x2p,u),x1,0),x2,0)];
C = [subs(subs(diff(x1,x1),x1,0),x2,0) subs(subs(diff(x2,x1),x1,0),x2,0)];
D = [0];

disp(' ')
disp('Matriz A = ')
disp(' ')
pretty(A);
disp(' ')
disp('Matriz B = ')
disp(' ')
pretty(B);
disp(' ')
disp('Matriz C = ')
disp(' ')
pretty(C);

m=2;
b=0.4;
l=1;
g=10;
delta=pi;
%---matrices
A =[0               , 1 ;
   (-g*cos(delta)/l), (-b/(l^2*m)) ]
B =[0 ; (1/(l^2*m))]

Ai =[0               , 1 ;
   (-g/l), (-b/(l^2*m)) ]
Bi =[0 ; (1/(l^2*m))]


eig(A)
eig(Ai)
%%
%  m=1.5;
% m=1.8;
% m=2;
m=2.2;
% m=2.8;
% m=3.3;
% m=3.8;
% m=4.4;
% m=4.9;
% m=5.5;
% m=6.6;

b=0.4;
delta=180; %en grados
l=1;
G=10;


% [A,B,C,D]=linmod('pendulo_mod_tarea',delta*pi/180)
% disp('Autovalores')
% eig(A)
% disp('Rango')
% rank(ctrb(A,B))
% Aa=[[A;C] zeros(3,1)]
% Ba=[B;0]
% disp('Autovalores')
% eig(Aa)
% disp('Rango')
% rank(ctrb(Aa,Ba))
% 
% 
% p= -2 %dato
% K=acker(Aa,Ba,[p p p])
% k1=K(1)
% k2=K(2)
% k3=K(3)
% disp('polos lazo cerrado')
% eig(Aa-Ba*K) % polos lazo cerrado
% disp('tiempo de respuesta calculado')
% tscalc=7.5/(-p) % tiempo de respuesta calculado

sim('pendulo_pid_tarea')
figure(1)
subplot(2,2,1);plot(tout,yout)
grid on,hold on, title('Salida')
subplot(2,2,2);, plot(yout,velocidad) %plano de fase
grid on,hold on, title('Plano de fases')
subplot(2,2,3);, plot(tout,torque) % torque total
grid on,hold on, title('Torque')
subplot(2,2,4);, plot(tout,-accint) % acción integral
grid on,hold on, title('Accion integral')
legend ('m=3.8','m=4.4','m=4.9','m=5.5','m=6.6');

ymax=max(yout) % máximo valor de salida
S=(ymax-delta)/delta*100 % sobrepaso en %
erel=(delta-yout)/delta; %error relativo
efinal=erel(end); % error final, debe ser cero
ind=find(abs(erel)>.02); % índice elementos con error relativo absoluto menor a 2%
tss=tout(ind(end)) % tiempo de establecimiento (ultimo valor del vector)
yte=yout(ind(end)) ;% salida al tiempo ts
uf=torque(end) ;% torque final
Intf=-accint(end) % acción integral final

