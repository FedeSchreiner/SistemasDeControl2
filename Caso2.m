%% punto 1 y 2

clear all,clc;
X=-[0 ; 0 ; 0 ];
%X=-[0 ; 0 ; 0 ];
ii=0;
t_etapa=1e-7;   %tiempo de integración
wRef=2;
tF=.15;        %tiempo de simulación
Ts=t_etapa;
e=zeros(tF/t_etapa,1);
u=12;
input=zeros(round(tF/t_etapa),1);
x1=0;
x2=0;
x3=0;
t1=0;

color_= 'k';
for t=0:t_etapa:tF
    ii=ii+1;
    X=modmotor(t_etapa, X, u); 
    x1(ii)=X(1); %Omega
    x2(ii)=X(2); %wp
    x3(ii)=X(3); %ia
    
    acc(ii)=u;
end

t=0:t_etapa:tF;

subplot(3,1,1);hold on;
plot(t,x1,color_);title('\omega_r');
subplot(3,1,2);hold on;
plot(t,x2,color_);title('i_a');
subplot(3,1,3);hold on;
% plot(t,x3,color_);title('Angulo \Theta');
% subplot(4,1,4);hold on;
plot(t,acc,color_);title('v_a');
xlabel('Tiempo [S]');



%% punto 3
%chema motor.
clear all;clc;close all;
%Sacamos los datos del excel
archivo = 'Curvas_Medidas_Motor1.xlsx';
hoja = 'Hoja1';

rango1= 'A3658:A6326';
rango2= 'B3658:B6326';
rango3= 'C3658:C6326';

t0=xlsread(archivo,hoja,rango1)-0.02;
omega=xlsread(archivo,hoja,rango2);
Ia=xlsread(archivo,hoja,rango3);

%%%aca comienza el algoritmo de chema
opt = stepDataOptions;
opt.StepAmplitude = 12;
t_inic=t0(1);

[val, lugar] =min(abs(t_inic-t0));
y_t1=omega(lugar);
t_t1=t0(lugar);

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

T1=T1_ang;
T2=T2_ang;
T3=T3_ang;

T3_ang=sum(T3/length(T3));
T2_ang=sum(T2/length(T2));
T1_ang=sum(T1/length(T1));

sys_G_ang=tf(K*[T3_ang 1],conv([T1_ang 1],[T2_ang 1]))
[y,t2] = step(sys_G_ang,opt);

%valores que me dio comparando con la formula.
Laa=13.93e-3;     %inductancia armadura
Ra=99;        %resistencia armadura
Ki=16.52;     %constantes del motor
Km=0.0605;     %contantes de morot
J=5.8441e-6;         %inercia motor
num=[Ki];
den=[Laa*J Ra*J Ki*Km];
G = tf(num/den(3),den/den(3))
[y1,t1] = step(G,opt);
figure(1);
step(sys_G_ang,opt);hold on;
step(G,opt),hold on;

%% -----------------------------
clear all;close all;clc
X=-[0 ; 0 ; 0;0];
ii=0;
t_etapa=1e-7;   %tiempo de integración
wRef=2;
tF=0.3;        %tiempo de simulación
Ts=t_etapa;
e=zeros(tF/t_etapa,1);
u=12;
input=zeros(round(tF/t_etapa),1);
delta_Tl=1e-7;
x1=0;
x2=0;
x3=0;
x4=0;
Tl=2e-5*0;
for t=0:t_etapa:tF
    ii=ii+1;
     if (ii==1e5)
         Tl=2e-5;
      end
    X=modmotor(t_etapa, X, u,Tl); %agregamos una varialbe de estado torqu
    x1(ii)=X(1); %Omega
    x2(ii)=X(2); %ia
    x3(ii)=X(3); %wp
    x4(ii)=X(4);
    acc(ii)=u;
end
t=0:t_etapa:tF;

% archivo = 'Curvas_Medidas_Motor.xlsx';
% hoja = 'Hoja1';
% rango1= 'A1:A22991';
% rango2= 'B1:B22991';
% rango3= 'C1:C22991';
% t0=xlsread(archivo,hoja,rango1);
% omega1=xlsread(archivo,hoja,rango3);
% Ia1=xlsread(archivo,hoja,rango2);
% subplot(2,1,1);hold on;
% plot(t,x1*0.11,color_);
% plot(t0,Ia1);
% title('\omega_r');
% subplot(2,1,2);hold on;
% plot(t,x2*0.31,color_);
% plot(t0,omega1);
% title('i_a');
% xlabel('Tiempo [S]');
% legend('modelo del chema','modelo del excel');

%% Punto 4

clear all;clc
X=-[0 ; 0 ; 0; 0 ];
ii=0;
t_etapa=1e-7;   %tiempo de integracion
wRef=2;
tF=0.5;        %tiempo de sumulacion

Kp=1e-1;
Ki=1e-3;
Kd=0;

color_='k';
Ts=t_etapa;
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;
e=zeros(tF/t_etapa,1);
u=0;
input=zeros(round(tF/t_etapa),1);
delta_Tl=1e-7;
Tl=0
for t=0:t_etapa:tF
    ii=ii+1;
    k=ii+2;
    if (ii==1e5)
          Tl=2e-5;
          end
    X=modmotor(t_etapa, X, u,Tl); %agregamos una variable de estado torque
    e(k)=wRef-X(4); %ERROR
    u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID
    x1(ii)=X(1); %Omega
    x2(ii)=X(2); %ia
    x3(ii)=X(3); %wp
    x4(ii)=X(4); %theta
    acc(ii)=u;
end
t=0:t_etapa:tF;
figure(1);
subplot(4,1,1);hold on;
plot(t,x4);title('\theta_t');
subplot(4,1,2);hold on;
plot(t,x1);title('\omega_t');
xlabel('Tiempo [Seg.]');
subplot(4,1,3);hold on;
plot(t,x3);title('Corriente');
subplot(4,1,4);hold on;
plot(t,acc);title('accion de control, u_t');
xlabel('Tiempo [Seg.]');

