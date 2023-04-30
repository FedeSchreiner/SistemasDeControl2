clear all;
close all;
clc;

%------constantes del motor---------
Laa=366e-6;
J=5e-9;         %inercia motor
Ra=55.6;        %resistencia armadura
Bi=0;            %amortiguamiento
Ki=6.49e-3;     %constantes del motor
Km=6.53e-3;     %contantes de motor

%------Matrices--------
Ap = [-Ra/Laa, -Km/Laa, 0;
     Ki/J,    -Bi/J   , 0;
     0,        1     , 0];
% Bc = [1/Laa, 0;
%      0, -1/J;
%      0,   0];
 Bp = [1/Laa;
     0 ;
     0   ];
Cp = [0,0,1];
Dp = [0];
% %D = [0, 0];

%-----Función de transferencia-----
[num,den] = ss2tf(Ap,Bp,Cp,Dp);
G = tf(num,den)

polos = eig(Ap)

%---tiempo muestreo
Ts=polos(2)

%-----Matrices discretas
sys_c=ss(Ap,Bp,Cp,Dp);
sys_d=c2d(sys_c,Ts,'zoh');

%-----Funcion de traferencia discreta-----
[num, den]=ss2tf(sys_d.A,sys_d.B,sys_d.C,sys_d.D,1);
Gd=tf(num,den);

%-----Función de transferencia discreta-----
[num, den]=ss2tf(sys_d.A,sys_d.B,sys_d.C,sys_d.D,1);
Gd=tf(num,den);
A = sys_d.A;
B = sys_d.B;
C = sys_d.C;
D = sys_d.D;

