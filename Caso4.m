%% Punto 1
clc;clear all;close all
clc;clear all;close all
m=0.1;
Fricc=0.1;
long=0.6;
g=9.8;
M=0.5;
h=1e-4;
tiempo=(250/h);
p_pp=0;
tita_pp=0;
t=0:h:tiempo*h;
omega=0:h:tiempo*h;
alfa=0:h:tiempo*h;
p=0:h:tiempo*h;
p_p=0:h:tiempo*h;
u=linspace(0,0,tiempo+1);
%Condiciones iniciales
alfa(1)=-0.01;
color='b';
p(1)=0;
p_p(1)=0;
u(1)=0;
p(1)=0;
i=1;
X0=[0 0 alfa(1) 0]';
x=[0 0 pi 0]';
for i=1:tiempo
    %Variables del sistema no lineaL
    u(i)=0;
    %Sistema no lineal
    p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))- Fricc*p_p(i));
    tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
    p_p(i+1)=p_p(i)+h*p_pp;
    p(i+1)=p(i)+h*p_p(i);
    omega(i+1)=omega(i)+h*tita_pp;
    alfa(i+1)=alfa(i)+h*omega(i);
end
figure(1);hold on;
subplot(4,1,1);plot(t,omega,color);grid on; title('Velocidad 聲gulo');hold on;%plot(t,omegal,'k');
subplot(4,1,2);plot(t,alfa,color);hold on;
% plot(t,pi*ones(size(t)),'k');
%plot(t,alfal,'k');
grid on;title('聲gulo');hold on;
subplot(4,1,3); plot(t,p,color);grid on;title('Posici鏮 carro');hold on;%plot(t,pl,'k');
subplot(4,1,4);plot(t,p_p,color);grid on;title('Velocidad carro');hold on;%plot(t,p_pl,'k');


%% punto 2
% lo mismo que en el punto anterior, pero con m=0.2 y M=1
clc;clear all;close all
m=0.1;
m1=m*2;
Fricc=0.1;
long=0.6;
g=9.8;
M=0.5;
M1=0.5*2;
h=1e-4;
tiempo=(6/h);
p_pp=0;
tita_pp=0;
p_pp1=0;
tita1_pp=0;
p_pp2=0;
tita2_pp=0;
t=0:h:(tiempo*h);
omega=0:h:tiempo*h;
alfa=0:h:tiempo*h;
p=0:h:tiempo*h;
p_p=0:h:tiempo*h;
omega1=0:h:tiempo*h;
alfa1=0:h:tiempo*h;
p1=0:h:tiempo*h;
p1_p=0:h:tiempo*h;
omega2=0:h:tiempo*h;
alfa2=0:h:tiempo*h;
p2=0:h:tiempo*h;
p2_p=0:h:tiempo*h;
u=linspace(0,0,tiempo+1);
%Condiciones iniciales
alfa(1)=-0.01;
alfa1(1)=-0.01;
alfa2(1)=-0.01;
% alfa(1)=3.01;
% alfa1(1)=3.01;
% alfa2(1)=3.01;
color='b';
p(1)=0;
p_p(1)=0;
u(1)=0;
p(1)=0;
i=1;
p1_p(1)=0;
p1(1)=0;
omega1(1)=0;
p2_p(1)=0;
p2(1)=0;
omega2(1)=0;
%Versi鏮 linealizada en el equilibrio estable. Sontag Pp 104.
X0=[0 0 alfa(1) 0]';
x=[0 0 pi 0]';
for i=1:tiempo
    %Variables del sistema no lineal
    %estado=[p(i); p_p(i); alfa(i); omega(i)];
    u(i)=0;
    %Sistema no lineal
    p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))- Fricc*p_p(i));
    tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
    p_p(i+1)=p_p(i)+h*p_pp;
    p(i+1)=p(i)+h*p_p(i);
    omega(i+1)=omega(i)+h*tita_pp;
    alfa(i+1)=alfa(i)+h*omega(i);
    p_pp1=(1/(M+m1))*(u(i)-m1*long*tita1_pp*cos(alfa1(i))+m1*long*omega1(i)^2*sin(alfa1(i))- Fricc*p1_p(i));
    tita1_pp=(1/long)*(g*sin(alfa1(i))-p_pp1*cos(alfa1(i)));
    p1_p(i+1)=p1_p(i)+h*p_pp1;
    p1(i+1)=p1(i)+h*p1_p(i);
    omega1(i+1)=omega1(i)+h*tita1_pp;
    alfa1(i+1)=alfa1(i)+h*omega1(i);
    p_pp2=(1/(M1+m))*(u(i)-m*long*tita2_pp*cos(alfa2(i))+m*long*omega2(i)^2*sin(alfa2(i))- Fricc*p2_p(i));
    tita2_pp=(1/long)*(g*sin(alfa2(i))-p_pp2*cos(alfa2(i)));
    p2_p(i+1)=p2_p(i)+h*p_pp2;
    p2(i+1)=p2(i)+h*p2_p(i);
    omega2(i+1)=omega2(i)+h*tita2_pp;
    alfa2(i+1)=alfa2(i)+h*omega2(i);
end
figure(1);hold on;
subplot(4,1,1);plot(t,omega,color);grid on; title('Velocidad 聲gulo');hold on;
plot(t,omega1,'k');hold on;
plot(t,omega2,'r');hold on;
subplot(4,1,2);plot(t,alfa,color);hold on;
plot(t,alfa1,'k');grid on
plot(t,alfa2,'r');grid on
grid on;title('聲gulo');hold on;
subplot(4,1,3); plot(t,p,color);grid on;title('Posici鏮 carro');hold on;
plot(t,p1,'k');
plot(t,p2,'r');
subplot(4,1,4);plot(t,p_p,color);grid on;title('Velocidad carro');hold on;
plot(t,p1_p,'k');
plot(t,p2_p,'r');
legend('Original','Masa m doble', 'Masa M doble');


%% punto 3
clear all;%13:44 16/1/2018
syms fi fi_p fi_pp p p_p p_pp M m u long Fricc g;
disp('Para el equilibrio inestable')
ang_inic=0;
p_pp=(1/(M+m))*(u-m*long*fi_pp+m*long*fi_p^2*fi-Fricc*p_p); %Peque隳s angulos
% fi_pp=(1/long)*(g*(fi)-p_pp); %Peque隳s angulos para fi~0, sin(fi)~fi, cos(fi)~1
fi_pp=solve(fi_pp==(1/long)*(g*fi-p_pp),fi_pp);
%disp('fi_pp='); pretty(simplify(fi_pp));
p_pp=subs(p_pp,'fi_pp',fi_pp);
Mat_A=[[0 1 0 0];
[subs(subs(subs(subs(diff(p_pp, p), p,0),p_p,0),fi,ang_inic),fi_p,0), ...
subs(subs(subs(subs(diff(p_pp, p_p), p,0),p_p,0),fi,ang_inic),fi_p,0), ...
subs(subs(subs(subs(diff(p_pp, fi), p,0),p_p,0),fi,ang_inic),fi_p,0), ...
subs(subs(subs(subs(diff(p_pp, fi_p), p,0),p_p,0),fi,ang_inic),fi_p,0)];
[0 0 0 1];
[subs(subs(subs(subs(diff(fi_pp, p), p,0),p_p,0),fi,ang_inic),fi_p,0),...
subs(subs(subs(subs(diff(fi_pp, p_p), p,0),p_p,0),fi,ang_inic),fi_p,0),...
subs(subs(subs(subs(diff(fi_pp, fi), p,0),p_p,0),fi,ang_inic),fi_p,0),...
subs(subs(subs(subs(diff(fi_pp, fi_p), p,0),p_p,0),fi,ang_inic),fi_p,0)]];
Mat_B=[0;
subs(subs(subs(subs(diff(p_pp, u), p,0),p_p,0),fi,ang_inic),fi_p,0);...
0;
subs(subs(subs(subs(diff(fi_pp, u), p,0),p_p,0),fi,ang_inic),fi_p,0)];
disp('Matriz A')
disp('  ')
pretty(simplify(Mat_A))
disp('Matriz B')
disp('  ')
pretty(simplify(Mat_B))


%% Punto 4
clc;clear all;close all
m=0.1;
m1=0.1;
Fricc=0.1;
long=1.2;
long1=1.2;
g=9.8;
M=0.5;
M1=0.5;
h=1e-4;
tiempo=(5/h);
p_pp=0;
tita_pp=0;
p_pp1=0;
tita1_pp=0;
p_pp2=0;
tita2_pp=0;
t=0:h:tiempo*h;
omega=0:h:tiempo*h;
alfa =0:h:tiempo*h;
p=0:h:tiempo*h;
p_p=0:h:tiempo*h;
omega1=0:h:tiempo*h;
alfal =0:h:tiempo*h;
pl = 0:h:tiempo*h;
p_pl=0:h:tiempo*h;
omega2=0:h:tiempo*h;
alfa2=0:h:tiempo*h;
p2=0:h:tiempo*h;
p2_p=0:h:tiempo*h;
u=linspace(0,0,tiempo+1);
%Condiciones iniciales
alfa(1)=0.8;
color='b';
p(1)=0;
p_p(1)=0;
u(1)=0;
p(1)=0;
i=1;

p_pl(1)=0;
pl(1)=0;

omega1(1)=0;
%Versi鏮 linealizada en el equilibrio estable. Sontag Pp 104.
%estado=[p(i); p_p(i); alfa(i); omega(i)]
Mat_A=[0 1 0 0;0 -Fricc/M -g*m/M 0;0 0 0 1; 0 Fricc/(M*long) g*(M+m)/(M*long) 0];
Mat_B=[0; 1/M; 0; -1/(M*long)];
X0=[0 0 0 0]';
x=[0 0 alfa(1) 0]';
for i=1:tiempo
%Variables del sistema no lineal
%estado=[p(i); p_p(i); alfa(i); omega(i)];
u(i)=0;
%Sistema no lineal
p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))- Fricc*p_p(i));
tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
p_p(i+1)=p_p(i)+h*p_pp;
p(i+1)=p(i)+h*p_p(i);
omega(i+1)=omega(i)+h*tita_pp;
alfa(i+1)=alfa(i)+h*omega(i);
 %Variables del sistema lineal
pl(i)=x(1);
p_pl(i)=x(2);
alfal(i)=x(3);
omegal(i+1)=x(4);
 %Sistema lineal
xp=Mat_A*(x-X0)+Mat_B*u(i);
x=x+h*xp;    
end
pl(i)=x(1);
p_pl(i)=x(2);
alfal(i)=x(3);
omegal(i)=x(4);
figure(1);hold on;
subplot(4,1,1);plot(t,omega,color);grid on; title('Velocidad 聲gulo');hold on;
plot(t,omegal,'r');hold on;
% plot(t,omega2,'r');hold on;
subplot(4,1,2);plot(t,alfa,color);hold on;
% plot(t,pi*ones(size(t)),'k');
plot(t,alfal,'r');grid on
% plot(t,alfa2,'r');grid on
grid on;title('聲gulo');hold on;
subplot(4,1,3); plot(t,p,color);grid on;title('Posici鏮 carro');hold on;
plot(t,pl,'r');
% plot(t,p2,'r');
subplot(4,1,4);plot(t,p_p,color);grid on;title('Velocidad carro');hold on;
plot(t,p_pl,'r');
% plot(t,p2_p,'r');
legend('no lineal','lineal');

%% punto 5

clear all;
syms fi fi_p fi_pp p p_p p_pp M m u long Fricc g;
disp('  ')
disp('Para el Pendulo en equilibrio estable')
disp('  ')
ang_inic=pi;
p_pp=(1/(M+m))*(u+m*long*fi_pp-m*long*fi_p^2*fi-Fricc*p_p);
%fi_pp=(1/long)*(-g*(fi)+p_pp); %Peque隳s angulos para fi~pi sin(fi)~-fi, cos(fi)=-1
fi_pp=solve(fi_pp==(1/long)*(-g*fi+p_pp),fi_pp);
disp('fi_pp='); pretty(simplify(fi_pp));
disp('  ')
p_pp=subs(p_pp,'fi_pp',fi_pp);
Mat_A=[[0 1 0 0];
[subs(subs(subs(subs(diff(p_pp, p), p,0),p_p,0),fi,ang_inic),fi_p,0), ...
subs(subs(subs(subs(diff(p_pp, p_p), p,0),p_p,0),fi,ang_inic),fi_p,0), ...
subs(subs(subs(subs(diff(p_pp, fi), p,0),p_p,0),fi,ang_inic),fi_p,0), ...
subs(subs(subs(subs(diff(p_pp, fi_p), p,0),p_p,0),fi,ang_inic),fi_p,0)]; ...
[0 0 0 1];...
[subs(subs(subs(subs(diff(fi_pp, p), p,0),p_p,0),fi,ang_inic),fi_p,0),...
subs(subs(subs(subs(diff(fi_pp, p_p), p,0),p_p,0),fi,ang_inic),fi_p,0),...
subs(subs(subs(subs(diff(fi_pp, fi), p,0),p_p,0),fi,ang_inic),fi_p,0),...
subs(subs(subs(subs(diff(fi_pp, fi_p), p,0),p_p,0),fi,ang_inic),fi_p,0)]];
Mat_B=[0;
subs(subs(subs(subs(diff(p_pp, u), p,0),p_p,0),fi,ang_inic),fi_p,0);...
0;
subs(subs(subs(subs(diff(fi_pp, u), p,0),p_p,0),fi,ang_inic),fi_p,0)];
disp('Matriz A')
disp('  ')
pretty(simplify(Mat_A))
disp('Matriz B')
disp('  ')
pretty(simplify(Mat_B))


%% punto 6

igual que el 4. cambia solo lo siguiente
alfa(1)=pi-0.8;



Mat_A=[0 1 0 0;0 -Fricc/M -g*m/M 0;0 0 0 1; 0 -Fricc/(M*long) -g*(M+m)/(M*long) 0];
Mat_B=[0; 1/M; 0; 1/(M*long)];
X0=[0 0 pi 0]';
x=[0 0 alfa(1) 0]';

%% c鏚igo extra para ver matriz no lineal.
clear all;clc;close all;
syms tita tita_p tita_pp p p_p p_pp M m u long Fricc g alfa;
disp('  ')
disp('Con X = [0 0 0.1 0] = No lineal')
disp('  ')
ang_inic=0;
p_pp=(1/(M+m))*(u-m*long*tita_pp*cos(tita)+m*long*(tita_p^2)*sin(tita)-Fricc*p_p);
tita_pp=(1/long)*(g*sin(tita)-p_pp*cos(tita));
Mat_A=[[0 1 0 0];
[subs(subs(subs(subs(diff(p_pp, p), p,0),p_p,0),fi,ang_inic),tita_p,0), ...
subs(subs(subs(subs(diff(p_pp, p_p), p,0),p_p,0),fi,ang_inic),tita_p,0), ...
subs(subs(subs(subs(diff(p_pp, tita), p,0),p_p,0),fi,ang_inic),tita_p,0), ...
subs(subs(subs(subs(diff(p_pp, tita_p), p,0),p_p,0),fi,ang_inic),tita_p,0)]; ...
[0 0 0 1];...
[subs(subs(subs(subs(diff(tita_pp, p), p,0),p_p,0),fi,ang_inic),tita_p,0),...
subs(subs(subs(subs(diff(tita_pp, p_p), p,0),p_p,0),fi,ang_inic),tita_p,0),...
subs(subs(subs(subs(diff(tita_pp, tita), p,0),p_p,0),fi,ang_inic),tita_p,0),...
subs(subs(subs(subs(diff(tita_pp, tita_p), p,0),p_p,0),fi,ang_inic),tita_p,0)]];
Mat_B=[0;
subs(subs(subs(subs(diff(p_pp, u), p,0),p_p,0),tita,ang_inic),tita_p,0);...
0;
subs(subs(subs(subs(diff(tita_pp, u), p,0),p_p,0),tita,ang_inic),tita_p,0)];
disp('  ')
disp('MATRIZ A')
disp('  ')
pretty(simplify(Mat_A))
disp('  ')
disp('MATRIZ B')
disp('  ')
pretty(simplify(Mat_B))


