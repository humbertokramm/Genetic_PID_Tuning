% Referência:
%http://ctms.engin.umich.edu/CTMS/index.php?example=Introduction&section=ControlPID

clc

%Parâmetros gerados no AG
Kp = 99.59349;
Ki = 98.81939;
Kd = 99.98095;

%Modelo do Controlador
C = Kp + Ki/s + Kd*s;
C = pid(Kp,Ki,Kd);

%Modelo da Planta
wn = 1;
alpha = 1;
csi = 0.25;
num = alpha*(wn^2);
den = s^2+(2*csi*wn)*s+(wn^2);

s = tf('s');
P = num/den
step(P)


C = pid(Kp,Ki,Kd);
T = feedback(C*P,1);

t = 0:0.01:2;
figure(1);
step(T,t)

%Tunador
opts = pidtuneOptions('CrossoverFrequency',32,'PhaseMargin',90);
[C, info] = pidtune(P, 'pid', opts);

%Resposta final
C
T = feedback(C*P,1);
t = 0:0.01:2;
figure(2);
step(T,t)
