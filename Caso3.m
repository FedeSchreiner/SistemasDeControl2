%%punto 1
clear all; clc;
syms alfa alfa_p p p_p p_pp h h_p u a b c omega;
disp('     ')
disp('Para el avion en equilibrio estable')
disp('     ')
%ecuaciones
alfa_p=a*(p-alfa);
p_pp=(-(omega^2))*(p-alfa-b*u);
h_p=c*alfa;
%matrices
Mat_A=[[subs(subs(subs(subs(diff(alfa_p, alfa), alfa,0),p,0),p_p,0),h,0), ...
 subs(subs(subs(subs(diff(alfa_p, p), alfa,0),p,0),p_p,0),h,0), ...
 subs(subs(subs(subs(diff(alfa_p, p_p), alfa,0),p,0),p_p,0),h,0), ...
 subs(subs(subs(subs(diff(alfa_p, h), alfa,0),p,0),p_p,0),h,0)]; ...
 [0 0 1 0];...
[subs(subs(subs(subs(diff(p_pp, alfa), alfa,0),p,0),p_p,0),h,0), ...
 subs(subs(subs(subs(diff(p_pp, p), alfa,0),p,0),p_p,0),h,0), ...
 subs(subs(subs(subs(diff(p_pp, p_p), alfa,0),p,0),p_p,0),h,0), ...
 subs(subs(subs(subs(diff(p_pp, h), alfa,0),p,0),p_p,0),h,0)]; ...
[subs(subs(subs(subs(diff(h_p, alfa), alfa,0),p,0),p_p,0),h,0), ...
 subs(subs(subs(subs(diff(h_p, p), alfa,0),p,0),p_p,0),h,0), ...
 subs(subs(subs(subs(diff(h_p, p_p), alfa,0),p,0),p_p,0),h,0), ...
 subs(subs(subs(subs(diff(h_p, h), alfa,0),p,0),p_p,0),h,0)]];
Mat_B=[0;
0;
subs(subs(subs(subs(diff(p_pp, u), alfa,0),p,0),p_p,0),h,0);...
0];
disp('MATRIZ A')
disp('  ')
pretty(simplify(Mat_A))
disp('  ')
disp('MATRIZ B')
disp('  ')
pretty(simplify(Mat_B))

%% Punto 2

clear all; clc ; close all;
syms p p_p p_pp alfa alfa_p ha ha_p a b c omega;
a=0.05;
b=5;
c=100;
omega=2;
color='b';
h=0.01;     %muestreo 10r-3
At=100;       %simulaci
tiempo=(At/h);   
alfa_p=0;
p_p=0;
p_pp=0;
ha_p=0;
t=0:h:tiempo*h;
alfa=0:h:tiempo*h;
p=0:h:tiempo*h;
p_p=0:h:tiempo*h;
ha=0:h:tiempo*h;
u=linspace(0,0,tiempo+1);
%condiciones iniciales
alfa(1)=0;
p(1)=0;
p_p(1)=0;
ha(1)=0;
u(1)=0;
i=0;
alfal(1)=0;
pl(1)=0;
p_pl(1)=0;
hal(1)=0;
%linealizada en eq. estable
Mat_A = [-a a 0 0; 0 0 1 0; omega^2 -omega^2 0 0; c 0 0 0];
Mat_B = [0;0;b*(omega^2);0];
x0=[0 0 0 0]';
x=[0 0 0 0]';
for i=1:tiempo
    if(i<=1000)
            u(i)= 0.2;
        end
        if(i>1000 && i<5000)
            u(i)= 0;
        end
         if(i>=5000 && i<6000)
            u(i)= -0.2;
         end
        if(i>=6000&& i<10000)
            u(i)=0;
        end
        %ecuaciones del sistema
        alfa_p=a*(p(i)-alfa(i));
        p_pp=(-omega^2)*(p(i)-alfa(i)-b*u(i));
        ha_p= c*alfa(i);
        alfa(i+1)=alfa(i)+h*alfa_p;
        p_p(i+1)=p_p(i)+h*p_pp;
        p(i+1)=p(i)+ h*p_p(i);
        ha(i+1)=ha(i)+h*ha_p;
        xp=Mat_A*(x-x0)+Mat_B*u(i);
            x=x+h*xp;
        %varibles y sistema lineal
            alfal(i+1)=x(1);
        pl(i+1)=x(2);
        p_pl(i+1)=x(3);
        hal(i+1)=x(4);
end
figure(1);hold on;
subplot(3,2,1);plot(t,alfal,color);grid on;title('Ángulo respecto al suelo \alpha');hold on;
subplot(3,2,2);plot(t,pl,color);hold on;
plot(t,pi*ones(size(t)),'k');grid on;title('Direccion avión \Phi');hold on;
subplot(3,2,3);plot(t,p_pl,color);hold on;
plot(t,pi*ones(size(t)),'k');grid on;title('Variacion de dirección');hold on;
subplot(3,2,4);plot(t,hal,color);grid on;title('Altura h');;hold on;
subplot(3,2,5);plot(t,u,color);grid on;title('u');;hold on;
color='k';
a=0.05; be=5; c=50 ; omega=2;   
%Paso, tiempo de simulacion y pasos;
At=10e-3; tF=20; pasos=tF/At; % aumentamos el timpo de simuación a 20
%Asignaciones
t=[]; u=linspace(0,0,pasos);
%Matriz de estados.
Mat_A=[-a a 0 0; 0 0 1 0; (omega^2) -(omega^2) 0 0 ; c 0 0 0];
%Matriz de entrada
Mat_B=[0; 0; (omega^2)*be; 0];
%Condiciones iniciales;
x=[0 0 0 0]'; X0=[0 0 0 0]'; u(1)=1; ii=1;
 
t=zeros(pasos,1);
alfal=zeros(pasos,1);
fil=zeros(pasos,1);
fi_pl=zeros(pasos,1);
hl=zeros(pasos,1);

while(ii<pasos+1)
    t(ii)=ii*At;
    u(ii)=1;
    %Variables del sistema lineal.
    alfal(ii)=x(1); fil(ii)=x(2);  fi_pl(ii)=x(3);  hl(ii)=x(4);
    %Sistema lineal
    xp=Mat_A*(x-X0)+Mat_B*u(ii);
    x=x+At*xp;
    ii=ii+1;
end
figure(1);
subplot(3,2,1);plot(t,alfal,color);grid on;title('Ángulo respecto al suelo \alpha');hold on;
subplot(3,2,2);plot(t,pl,color);hold on;
plot(t,pi*ones(size(t)),'k');grid on;title('Direccion avión \Phi');hold on;
subplot(3,2,3);plot(t,p_pl,color);hold on;
plot(t,pi*ones(size(t)),'k');grid on;title('Variacion de dirección');hold on;
subplot(3,2,4);plot(t,hal,color);grid on;title('Altura h');;hold on;
subplot(3,2,5);plot(t,u,color);grid on;title('u');;hold on;

%% Punto 3
%lo mismo que en punto anterior pero como c=50
