% Referência:
%http://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=ControlPID
clear
clc
s = tf('s');

%Parâmetros gerados no AG
Kp = 99.59349;
Ki = 98.81939;
Kd = 99.98095;

%Modelo do Controlador
C = Kp + Ki/s + Kd*s;
C = pid(Kp,Ki,Kd);

%Modelo da Planta 1
wn = 1;
alpha = 1;
csi = 0.25;
num = alpha*(wn^2);
den = s^2+(2*csi*wn)*s+(wn^2);

s = tf('s');
P = num/den

step(P)

hold on

C = pid(Kp,Ki,Kd);
T_ga = feedback(C*P,1);
    
t = 0:0.01:5;
step(T_ga,t)

% parametros do tuner: 
% Kp = 96.82
% Ki = 44.58
% Kd = 52.58

%Tunador
opts = pidtuneOptions('CrossoverFrequency',52.61,'PhaseMargin',90);
[C, info] = pidtune(P, 'pid', opts);

%Resposta final
C
T_matlab = feedback(C*P,1);
t = 0:0.01:5;
step(T_matlab,t)
