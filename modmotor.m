function [X]=modmotor(t_etapa,xant,accion)
Laa=366e-6;     %inductancia armadura
J=5e-9;         %inercia motor
Ra=55.6;        %resistencia armadura
B=0;            %amortiguamiento
Ki=6.49e-3;     %constantes del motor
Km=6.53e-3;     %contantes de motot

TL=2.12788e-5;
Va=accion;   %voltaje de entrada
h=1e-7;      
omega=xant(1);
ia=xant(2);
wp=xant(3);
    for ii=1:t_etapa/h
        wpp= (-wp*(Ra*J+Laa*B)-omega*(Ra*B+Ki*Km)+Va*Ki)/(J*Laa);
        iap= ((-Ra*ia)-(Km*omega)+Va)/Laa;
        wp= wp+h*wpp;
        wp= wp -((1/J)*TL);
        omega=omega+h*wp;
        ia=ia+iap*h;
    end
X=[omega,ia,wp];


%motor:
% function [X]=modmotor(t_etapa,xant,accion,Tl)
% Laa=366e-6;
% J=5e-9;         %inercia motor
% Ra=55.6;        %resistencia armadura
% B=0;            %amortiguamiento
% Ki=6.49e-3;     %constantes del motor
% Km=6.53e-3;     %contantes de morot
% TL=Tl;
% Va=accion;   %voltaje de entrada
% h=1e-7;      % como si fuera el tiempo de simulacion
% omega=xant(1);
% wp=xant(2);
% ia=xant(3);
% thete=xant(4);
%     for ii=1:t_etapa/h
%         wpp= (-wp*(Ra*J+Laa*B)-omega*(Ra*B+Ki*Km)+Va*Ki)/(J*Laa);
%         iap= ((-Ra*ia)-(Km*omega)+Va)/Laa;
%         wp= wp+h*wpp;
%         omega=omega+h*wp;
%         ia=ia+iap*h;
%           thetep = omega;
%         thete = thete + h*thetep;
%       end
% X=[omega,wp,ia,thete];
