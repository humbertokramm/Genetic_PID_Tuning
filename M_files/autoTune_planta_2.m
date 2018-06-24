% Referência:
%http://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=ControlPID
clear
clc
s = tf('s');

%Parâmetros gerados no AG

Kp = 99.94621 
Ki = 30.98894 
Kd = 99.90348 

%Modelo do Controlador
C = Kp + Ki/s + Kd*s;
C = pid(Kp,Ki,Kd);

%Modelo da Planta 2
wn = 0.6;
alpha = 1;
csi = 0.1;
num = alpha*(wn^2);
den = s^2+(2*csi*wn)*s+(wn^2);

s = tf('s');
P = num/den

step(P)

hold on

C = pid(Kp,Ki,Kd);
T_ga = feedback(C*P,1);

t = 0:0.01:30;
step(T_ga,t)


% parametros do tuner: 
% Kp = 97.6475
% ki = 27.0883
% Kd = 97.9995

%Tunador
opts = pidtuneOptions('CrossoverFrequency',31.61,'PhaseMargin',90);
[C, info] = pidtune(P, 'pid', opts);

%Resposta final
C
T_matlab = feedback(C*P,1);
t = 0:0.01:30;
step(T_matlab,t)
