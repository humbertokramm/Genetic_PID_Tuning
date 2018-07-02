% Parâmetros do sistema
wn = 1;
csi = 0.25;
alpha = 1;

% Parâmetros da resposta
tMax = 20;      % Tempo final de amostragem
tRes = 0.001;   % Resolução do tempo de amostragem

% Cria sistema de segunda ordem
sys = tf([alpha*wn], [1 2*csi*wn wn*wn]);

% Obtém a resposta do sistema ao degrau (entre 0s e 20s em passos de 1ms)
[x, y] = step(sys, 0:0.001:20);
step(sys, 0:tRes:tMax);

% Escreve o CSV com tempo e resposta do sistema
csvwrite('resposta.csv', [x(:), y(:)]);