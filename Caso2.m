%% punto 1 y 2

clear all,clc;
X=-[0 ; 0 ; 0 ];

ii=0;
t_etapa=1e-7;   %tiempo de integración
wRef=2;
% tF=.08;        %tiempo de simulación
 tF=5;
Ts=t_etapa;
e=zeros(tF/t_etapa,1);
u=12;
input=zeros(round(tF/t_etapa),1);
x1=0;
x2=0;
x3=0;
t1=0;

color_= 'g';
for t=0:t_etapa:tF
    ii=ii+1;
    X=modmotor(t_etapa, X, u); 
    x1(ii)=X(1); %Omega
    x2(ii)=X(2); %wp
    x3(ii)=X(3); %ia
    
    acc(ii)=u;
end

t=0:t_etapa:tF;

% subplot(3,1,1);hold on;
% plot(t,x1,color_);title('\omega_r');xlabel('Tiempo [S]');
% subplot(3,1,2);hold on;
% plot(t,x2,color_);title('i_a');xlabel('Tiempo [S]');
% subplot(3,1,3);hold on;
% % plot(t,x3,color_);title('Angulo \Theta');
% % subplot(4,1,4);hold on;
% plot(t,acc,color_);title('v_a');xlabel('Tiempo [S]');

plot(t,x2);title('i_a');xlabel('Tiempo [S]');ylabel('[A]');grid on;hold on;
 legend('Con TL=0','Con TL=2e-7','Con TL=2e-6','Con TL=2.12788e-5');


%% punto 3


%chema motor.
clear all;clc;close all
%Sacamos los datos del excel
archivo = 'Curvas_Medidas_Motor_2023.xlsx';
hoja = 'Hoja1';

% rango1= 'A3658:A6326';
% rango2= 'B3658:B6326';
% rango3= 'C3658:C6326';

% rango1= 'A101:A31054';
% rango2= 'B101:B31054';
% rango3= 'C101:C31054';
% rango4= 'D101:D31054';
% rango5= 'E101:E31054';

% sin torque
rango1= 'A101:A15306';
rango2= 'B101:B15306';
rango3= 'C101:C15306';
rango4= 'D101:D15306';
rango5= 'E101:E15306';
% 
% %con torque
% rango1= 'A15307:A31054';
% rango2= 'B15307:B31054';
% rango3= 'C15307:C31054';
% rango4= 'D15307:D31054';
% rango5= 'E15307:E31054';

% t0=xlsread(archivo,hoja,rango1)-0.1504;
t0=xlsread(archivo,hoja,rango1)-0.025;
omega=xlsread(archivo,hoja,rango2);
Ia=xlsread(archivo,hoja,rango3);
Va=xlsread(archivo,hoja,rango4);
Tll=xlsread(archivo,hoja,rango5);

%%%aca comienza el algoritmo de chema
opt = stepDataOptions;
opt.StepAmplitude = 12;
t_inic=t0(1669);

[val, lugar] =min(abs(t_inic-t0));
y_t1=omega(lugar);
t_t1=t0(lugar);

ii=0; 
ii=ii+1;

[val, lugar] =min(abs(2*t_inic-t0));
t_2t1=t0(lugar);
y_2t1=omega(lugar);

[val, lugar] =min(abs(3*t_inic-t0));
t_3t1=t0(lugar);
y_3t1=omega(lugar);

K=omega(end)/opt.StepAmplitude;

k1=(1/opt.StepAmplitude)*y_t1/K-1;
k2=(1/opt.StepAmplitude)*y_2t1/K-1;
k3=(1/opt.StepAmplitude)*y_3t1/K-1;

be=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1=(k1*k2+k3-sqrt(be))/(2*(k1^2+k2));
alfa2=(k1*k2+k3+sqrt(be))/(2*(k1^2+k2));

beta=(k1+alfa2)/(alfa1-alfa2);

T1_ang=-t_t1/log(alfa1);
T2_ang=-t_t1/log(alfa2);
T3_ang=beta*(T1_ang-T2_ang)+T1_ang;

T1(ii)=T1_ang;
T2(ii)=T2_ang;
T3(ii)=T3_ang;

T3_ang=sum(T3/length(T3));
T2_ang=sum(T2/length(T2));
T1_ang=sum(T1/length(T1));

sys_G_ang=tf(K*[T3_ang 1],conv([T1_ang 1],[T2_ang 1]))
[y,t2] = step(sys_G_ang,opt);

%valores que me dio comparando con la formula.
Laa=620e-6;     %inductancia armadura
Ra=99.24;        %resistencia armadura
Ki=16.52;     %constantes del motor
Km=0.0605;     %contantes de morot
J=4.1748e-6;         %inercia motor
num=[Ki];
den=[Laa*J Ra*J Ki*Km];
G = tf(num/den(3),den/den(3))
[y1,t1] = step(G,opt);

figure(1);
step(sys_G_ang,opt);hold on;

step(G,opt),hold on;
legend('modelo original','modelo obtenido');
figure(2);
plot(t0,omega);hold on;grid on;
subplot(4,1,1),plot(t0,omega);ylabel('Wr[rad s^-1]');grid on;hold on;
subplot(4,1,2),plot(t0,Ia);ylabel('Ia[A]');grid on;hold on;
subplot(4,1,3),plot(t0,Va);ylabel('Va [V]');grid on;hold on;
subplot(4,1,4),plot(t0,Tll);xlabel('Tiempo [S]');ylabel('TL[Nm]');grid on;hold on;

%% -----------------------------
clear all;close all;clc
X=-[0 ; 0 ; 0;0];
ii=0;
k=0;
t_etapa=1e-7;   %tiempo de integración

tF=0.6;        %tiempo de simulación

Ts=t_etapa;
e=zeros(tF/t_etapa,1);
u=12;
input=zeros(round(tF/t_etapa),1);
delta_Tl=1e-7;
x1=0;
x2=0;
x3=0;
x4=0;
Tl=7.5e-5*0;
for t=0:t_etapa:tF
    ii=ii+1;
    k=k+1;
     if (ii==1.2510e+06)
         u=-12;
     end
      if (ii==1.2538e+06)
         Tl=-2.74e-2;
      end
    X=modmotor(t_etapa, X, u,Tl); %agregamos una varialbe de estado torqu
    x1(ii)=X(1); %Omega
    x2(ii)=X(2); %wp
    x3(ii)=X(3); %ia
%     x4(ii)=X(4);
    acc(ii)=u;
end
t=0:t_etapa:tF;

archivo = 'Curvas_Medidas_Motor_2023.xlsx';
hoja = 'Hoja1';
rango1= 'A101:A31054';
rango2= 'B101:B31054';
rango3= 'C101:C31054';

color_='b';
t0=xlsread(archivo,hoja,rango1)-0.025;
omega1=xlsread(archivo,hoja,rango2);
Ia1=xlsread(archivo,hoja,rango3);

figure(1)
subplot(2,1,1);hold on;
plot(t,x1,color_);
plot(t0,omega1);
title('\omega_r');

subplot(2,1,2);hold on;
plot(t,x3,color_);
plot(t0,Ia1);
title('i_a');
xlabel('Tiempo [S]');
legend('modelo del chema','modelo del excel');

%% Punto 4

clear all;clc;
X=-[0 ; 0 ; 0; 0 ];
ii=0;
t_etapa=1e-6;   %tiempo de integracion
wRef=1;
tF=0.15;        %tiempo de sumulacion

Kp=23;
% Kp=30;
% Ki=10;
% Kd=0.000000001;
% Kp=5;
% Ki=10;
Ki=25;
Kd=0.000001;

color_='k';
Ts=t_etapa;
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;
e=zeros(tF/t_etapa,1);
u=0;
input=zeros(round(tF/t_etapa),1);
% Tl=-2.74e-2
Tl=1e-4;
for t=0:t_etapa:tF
    ii=ii+1;
    k=ii+2;
%     if (ii==1e5)
%           Tl=-2.74e-2;
%           end
    X=modmotor(t_etapa, X, u,Tl); %agregamos una variable de estado torque
    e(k)=wRef-X(4); %ERROR
    u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID
    x1(ii)=X(1); %Omega
    x2(ii)=X(2); %wp
    x3(ii)=X(3); %ia
    x4(ii)=X(4); %theta
    acc(ii)=u;
end
t=0:t_etapa:tF;
figure(1);
subplot(1,4,1);hold on;
plot(t,x4);title('\theta_t');

subplot(1,4,2);hold on;
plot(t,x1);title('\omega_t');
xlabel('Tiempo [Seg.]');

subplot(1,4,3);hold on;
plot(t,x3);title('Corriente');

subplot(1,4,4);hold on;
plot(t,acc);title('accion de control, u_t');
xlabel('Tiempo [Seg.]');

% legend('','');
