clear all; close all; clc;

R = 4.7e3;
L = 10e-6;
C = 100e-9;

Mat_A = [-R/L -1/L; 1/C 0];
Mat_B = [ 1/L ; 0];
Mat_C = [R 0];
Mat_D = 0;

[num,den] = ss2tf(Mat_A,Mat_B,Mat_C,Mat_D);
G=tf(num,den)
a = eig(G)
% At = log(0.95)/ a(1) * 0.9   %tiempo de muestreo
% T = log(0.05)/ a(2) * 2.4    %tiempo de simulación
At = 36.38e-11   %tiempo de muestreo
T = 3e-3   %tiempo de simulación

%condiciones iniciales;
Vc(1)=0;
Il(1)=0;
u(1)=12;
Kmax=round(T/At)
I1=zeros(1,Kmax);
Vc=zeros(1,Kmax);
u=linspace(0,0,Kmax);
i=1;
ii=0;
vin=12;
t=linspace(0,T,Kmax);
x=[Il(1) Vc(1)]';    
X0=[0 0]';

for i=1:Kmax
    ii=ii+At;
    
       if(ii>=1e-3)
       vin=vin*(-1);
       ii=0;
       end
       
       %variables del sistema
    u(i)=vin;
    xp=Mat_A*(x-X0)+Mat_B*u(i);
    x=x+xp*At;
    Il(i)= x(1);
    Vc(i)= x(2);
    
end

subplot(3,1,1);hold on;grid on;plot(t,Il,'b');title('IL,r');xlabel('Tiempo [S]');ylabel('[A]');
subplot(3,1,2);hold on;grid on;plot(t,Vc,'c');title('Vc,r');xlabel('Tiempo [S]');ylabel('[V]');
subplot(3,1,3);hold on;grid on;plot(t,u,'m');title('U,r');xlabel('Tiempo [S]');ylabel('[V]');


%% Punto 2
clear all; close all; clc;
R = 5.6e3; 
L = 10e-6;
C = 100e-9;

Mat_A = [-R/L -1/L; 1/C 0];
Mat_B = [ 1/L ; 0];
Mat_C = [R 0];
Mat_D = 0;

[num,den] = ss2tf(Mat_A,Mat_B,Mat_C,Mat_D);
G=tf(num,den)
a = eig(G);

% At = (log(0.95)/ a(1)) ; %tiempo de muestreo
% T = (log(0.05)/ a(2)) ;  %tiempo de simulación

At = 3e-10;
T=3e-3;
% At = 36.38e-11   %tiempo de muestreo
% T = 3e-3   %tiempo de simulación

%condiciones iniciales;
Vc(1)=0;
Il(1)=0;
u(1)=12;
Kmax=round(T/At);

I1=zeros(1,Kmax);
Vc=zeros(1,Kmax);
u=linspace(0,0,Kmax);
i=1;
ii=0;
vin=12;
t=linspace(0,T,Kmax);
x=[Il(1) Vc(1)]';    
X0=[0 0]';
for i=1:Kmax
    ii=ii+At;
       if(ii>=1e-3)
           vin=vin*(-1);
           ii=0;
       end
       %variables del sistema
    u(i)=vin;
    xp=Mat_A*(x-X0)+Mat_B*u(i);
    x=x+xp*At;
    Il(i)= x(1);
    Vc(i)= x(2);
       end
subplot(3,1,1);hold on;grid on;plot(t,Il,'m');title('IL,r');xlabel('Tiempo [S]');ylabel('[A]');
subplot(3,1,2);hold on;grid on;plot(t,Vc,'b');title('Vc,r');xlabel('Tiempo [S]');ylabel('[V]');
subplot(3,1,3);hold on;grid on;plot(t,u,'c');title('U,r');xlabel('Tiempo [S]');ylabel('[V]');

%% PUNTO 3
% el modelo dinamico a traves del metodo de chena
clear all;clc;close all
%Sacamos los datos del excel
archivo = 'Curvas_Medidas_RLC.xls';
hoja = 'Hoja1';

%solo donde hace el step
rango1= 'A102:A501';
rango2= 'B102:B501';
rango3= 'C102:C501';

t0 = xlsread(archivo,hoja,rango1)-0.01; % le saco el retardo.
i  = xlsread(archivo,hoja,rango2);
Vc = xlsread(archivo,hoja,rango3);

%%%aca comienza el algoritmo de chema
opt = stepDataOptions;
opt.StepAmplitude = 12;
t_inic=t0(20); 

ii=0; % for t_inic=10:15
ii=ii+1;

[val, lugar] =min(abs(t_inic-t0));
y_t1=Vc(lugar);
t_t1=t0(lugar);


[val, lugar] =min(abs(2*t_inic-t0));
t_2t1=t0(lugar);
y_2t1=Vc(lugar);

[val, lugar] =min(abs(3*t_inic-t0));
t_3t1=t0(lugar);
y_3t1=Vc(lugar);


K=Vc(end)/opt.StepAmplitude;

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
[y,t0] = step(sys_G_ang,opt);


%planteamos un R que a nosotros nos parezca
R = 100
%tenemos que abajo para despejar
[numerador denominador] = tfdata(sys_G_ang, 'v');
denominador = denominador / denominador(1);
%ahora que tenemos normalizado el denominador, lo utilizamos
%para calcular L y C
L = R / denominador(2)
C = 1 /(L * denominador(3))
%normalizamos el numerador sin la componente de laplace
numerador = numerador / numerador(3);
%lo colocamos como la estructura
numerador(3) = numerador(3) * 1/(L*C);
G = tf(numerador,[1 (R/L) 1/(L*C)]);
[y1,t1] = step(G,opt);
figure(1);
step(sys_G_ang,opt);hold on;
step(G,opt),hold on;
legend('Respuesta al escalón Excel','Respuesta al escalón Modelo');



%% Punto 4
close all,clear all,clc;
%emplear la serie de corriente desde 0.05s en adelante para validar el
%resultado.

% Hay que hacer el escalon invertido y ver como nuestra eleccion de
% resistencia, inducto y capacitor nos da muy parecido a el excel.
% para esto, hay que graficar la corriente del excel.

%con los valores de RLC y la funcion de tranferencia de Voltaje de
%entrada,y salida corriente.

% R= 100;
% L= 0.045205424253276;
% C= 4.830605770781058e-05;

R= 110;
L= 0.0497;
C= 4.3915e-05;

Mat_A = [-R/L -1/L; 1/C 0];
Mat_B = [ 1/L ; 0];
Mat_C = [R 0];
Mat_D = 0;

[num,den] = ss2tf(Mat_A,Mat_B,Mat_C,Mat_D);
G=tf(num,den)
a = eig(G);

% At = (log(0.95)/ a(1)) * 1.5 ; %tiempo de muestreo
% T = (log(0.05)/ a(2)) *0.01 ;  %tiempo de simulación

At = 1e-5;
T=5e-2;

Vc(1)=0;
I(1)=0;
u(1)=12;
Kmax=round(T/At);
I=zeros(1,Kmax);
Vc=zeros(1,Kmax);
u=linspace(0,0,Kmax);
i=1;

vin=12;
t=linspace(0,T,Kmax);
x=[I(1) Vc(1)]';    
X0=[0 0]';;

for i=1:Kmax
    
        u(i)=vin;
        xp=Mat_A*(x-X0)+Mat_B*u(i);
        x=x+xp*At;
        I(i)= x(1);
        Vc(i)= x(2);
    
        vin=-12;
end


%----
archivo = 'Curvas_Medidas_RLC.xls';
hoja = 'Hoja1';
%solo donde hace el step
rango1= 'A501:A1001';
rango2= 'B501:B1001';

t0 = xlsread(archivo,hoja,rango1)-0.05; % lo ubico en el punto
im  = xlsread(archivo,hoja,rango2);

figure(1);
plot(t0,im,'k');hold on;
plot(t,I,'b');hold on;
legend('Respuesta al escalón Excel','Respuesta al escalón Modelo');

