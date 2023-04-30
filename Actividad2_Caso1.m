%%codigo Item [1] actividad 2
%solo 3 variables de estado. 

clear all;close all;clc

Laa=366e-6;
J=5e-9;         %inercia motor
Ra=55.6;        %resistencia armadura
B=0;            %amortiguamiento
Ki=6.49e-3;     %constantes del motor
Km=6.53e-3;     %contantes de morot

ii=0;
t_etapa=1e-7;  
tF=1;        
Kmax = tF/t_etapa;

t=linspace(0,t_etapa,Kmax);
u=linspace(0,t_etapa,Kmax);

%condiciones iniciales
omega(1)= 0;
ia(1)= 0;
thete(1)=0;
wp(1)= 0;
u(1)=0; 

TL=0;
TLpi=1.15e-3;


X=[omega(1) ia(1) thete(1) ]';
Xop=[0 0 0]';

% tp=0,05; %tiempo perturbacion

for i=0:Kmax/10
    ii=ii+1;
    
    %      if (ii1e5)
    % %             Tl=7.5e-2;
    %          Tl=2e-5;
    %       end


    %modelo en ecuaciones
    wpp= (-wp*(Ra*J+Laa*B)-omega*(Ra*B+Ki*Km)+u*Ki)/(J*Laa);
    iap= ((-Ra*ia)-(Km*omega)+u)/Laa;

    wp = wp+t_etapa*wpp;
    wp = wp -((1/J)*TL);
    omega=omega+t_etapa*wp;
    ia=ia+iap*t_etapa;
    thete = thete + t_etapa*omega;
end

% t1(i+1)=Tl;
% T1=T1+delta_T1;


figure(1);
subplot(4,1,1);hold on;
plot(t,x3);title('\theta_t');
subplot(4,1,2);hold on;
plot(t,x1);title('\omega_t');
xlabel('Tiempo [Seg.]');
subplot(4,1,3);hold on;
plot(t,x2);title('Corriente');
subplot(4,1,4);hold on;
plot(t,u);title('accion de control, u_t');
xlabel('Tiempo [Seg.]');

%% vamos a copiar el ejemplo del pendulo y llevarlo a lo del motor
clc;clear all;

%datos
m=.1;
Fricc=0.1;
long=0.6;
g=9.8;
M=.5;

%constantes de tiempo
h=0.001;
tiempo=(50/h);

p_pp=0;
tita_pp=0;

%Condiciones iniciales
alfa(1)=.1; color='r';
% alfa(1)=.2; color='g';
alfa(1)=.9; color='b';

omega(1)=0;
p_p(1)=0;
u(1)=0;
p(1)=0;
i=1;
indice=0;

% % estado=[p(i); p_p(i); alfa(i); omega(i)]
% Mat_A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(long*M) g*(m+M)/(long*M) 0]
% Mat_B=[0; 1/M; 0; -1/(long*M)]
% Mat_C=[1 0 0 0]; %La salida es posición y ángulo
% Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B Mat_A^3*Mat_B ];%Matriz Controlabilidad
% %Cálculo del controlador por asignación de polos
% auto_val=eig(Mat_A);
% c_ai=conv(conv(conv([1 -auto_val(1)],[1 -auto_val(2)]),[1 -auto_val(3)]),[1 -
% auto_val(4)]);
% Mat_W=[c_ai(4) c_ai(3) c_ai(2) 1;c_ai(3) c_ai(2) 1 0;c_ai(2) 1 0 0;1 0 0 0];
% Mat_T=Mat_M*Mat_W;
% A_controlable=inv(Mat_T)*Mat_A*Mat_T %Verificación de que T esté bien
% %CONTROLADOR Ubicación de los polos de lazo cerrado en mui :
% mui(1)=-20;mui(2)=-5; mui(3)=-.05 + 0.5i;mui(4)=conj(mui(3));
% alfa_i=conv(conv(conv([1 -mui(3)],[1 -mui(4)]),[1 -mui(2)]),[1 -mui(1)]);
% K=fliplr(alfa_i(2:5)-c_ai(2:5))*inv(Mat_T);
% eig(Mat_A-Mat_B*K)
% 
% Mat_A_O=Mat_A';
% Mat_B_O=Mat_C';
% Mat_M_Dual=[Mat_B_O Mat_A_O*Mat_B_O Mat_A_O^2*Mat_B_O Mat_A_O^3*Mat_B_O];%Matriz
% %Controlabilidad
% alfaO_i=alfa_i;
% % Ubicacion del Observador
% % Algunas veces más rápido que el controlador
% mui_o=real(mui)*20;
% alfaO_i=conv(conv(conv([1 -mui_o(3)],[1 -mui_o(4)]),[1 -mui_o(2)]),[1 -mui_o(1)]);
% Mat_T_O=Mat_M_Dual*Mat_W;
% Ko=(fliplr(alfaO_i(2:end)-c_ai(2:end))*inv(Mat_T_O))';
% eig(Mat_A_O'-Ko*Mat_C) %Verifico que todos los polos estén en el semiplano izquierdo
% x_hat=[0;0;0;0]; %Inicializo el Observador
while(i<(tiempo+1))
    
    estado=[p(i); p_p(i); alfa(i); omega(i)]; 
    
%     % u(i)=-K*estado; color='*-b'; %Sin Observador
%     u(i)=-K*x_hat; color='.-r'; %Con Observador
%     
    
    %ecuaciones
    p_pp=(1/(M+m))*(u(i)-m*long*tita_pp*cos(alfa(i))+m*long*omega(i)^2*sin(alfa(i))- Fricc*p_p(i));
    tita_pp=(1/long)*(g*sin(alfa(i))-p_pp*cos(alfa(i)));
   
    p_p(i+1)=p_p(i)+h*p_pp;
    p(i+1)=p(i)+h*p_p(i);
    omega(i+1)=omega(i)+h*tita_pp;
    alfa(i+1)=alfa(i)+h*omega(i);
    
%     y_sal(i)=Mat_C*estado;
%     %________OBSERVADOR__________
%     y_sal_O(i)=Mat_C*x_hat;
%     y_sal(i)=Mat_C*estado;
%     x_hatp=Mat_A*x_hat+Mat_B*u(i)+Ko*(y_sal(i)-y_sal_O(i));
%     x_hat=x_hat+h*x_hatp;
%     i=i+1;
end
figure(1);hold on; t=1:i;t=t*h;
subplot(3,2,1);plot(t,omega,color);grid on; title('Velocidad ángulo');hold on;
legend('Sin Observador','Con Observador');legend('boxoff');
subplot(3,2,2);plot(t,alfa,color);grid on;title('Ángulo');hold on;
subplot(3,2,3); plot(t,p,color);grid on;title('Posición carro');hold on;
subplot(3,2,4);plot(t,p_p,color);grid on;title('Velocidad carro');hold on;
subplot(3,1,3);plot(t(1:end-1),u,color);grid on;title('Acción de control');xlabel('Tiempo en Seg.');hold on;
figure(2);hold on;
subplot(2,2,1);plot(alfa,omega,color);grid on;xlabel('Ángulo');ylabel('Velocidad angular');hold on;
subplot(2,2,2);plot(p,p_p,color);grid on;xlabel('Posicion carro');ylabel('Velocidad carro');hold on;
legend('Sin Observador','Con Observador');legend('boxoff');
