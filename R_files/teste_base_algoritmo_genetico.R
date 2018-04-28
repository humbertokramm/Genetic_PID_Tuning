####### Configuracoes do programa #######

# Modelo do sistema de segunda ordem (as variaveis importantes)
alpha = 4
wn = 1
csi = 0.5

# Configuracoes do algoritmo genetico
geracoes = 100
populacao = 50

Kp_min = 0
Kp_max = 100

Ki_min = 0
Ki_max = 100

Kd_min = 0
Kd_max = 100

tx_amostragem = 5000

#########################################

# Modelo do sistema de segunda ordem
sistema = function(s) (alpha*(wn^2))/(s^2 + 2*csi*wn*s + wn^2)

# Biblioteca que contem invlap() - transformada inversa de laplace
library(pracma)
library(genalg)

# Funcao que avalia o cromossomo
# cromossomo = (P, I, D)
fitness = function(cromossomo) {
  # Separa o cromossomo nas variaveis de ajuste do PID
  Kp = cromossomo[1]
  Ki = cromossomo[2]
  Kd = cromossomo[3]
  
  # Modelo do controlador PID
  pid = function(s) (Kp + (Ki/s) + Kd*s)
  
  # PID e sistema com realimentacao recebendo entrada de um degrau
  malha_fechada = function(s) (pid(s)*sistema(s))/(1 + pid(s)*sistema(s))

  # Resposta do sistema realimentado ao degrau unitario no dominio da frequencia
  resposta_degrau_unitario_f = function(s) (1/s)*malha_fechada(s)
  
  # Resposta do sistema com controlador ao degrau unitário no dominio do tempo
  resposta_do_sistema_t = invlap(resposta_degrau_unitario_f, 0, 2*pi, tx_amostragem)
  
  # A funcao de fitness retorna a soma do erro quadratico em relacao a uma exponencial (resposta utopica)
  erro = 0
  for (i in 1:(tx_amostragem-1)){
    y_desejado = (1-exp(-10*resposta_do_sistema_t$x[i]))
    erro = erro + (y_desejado-resposta_do_sistema_t$y[i])^2
  }
  return(erro)
}

# Algoritmo genetico (https://www.rdocumentation.org/packages/genalg/versions/0.2.0/topics/rbga)
GAmodel = rbga(
  # Valores m�xnimos de cada componente do cromossomo
  stringMin=c(Kp_min, Ki_min, Kd_min),
  # Valores maximos de cada componente do cromossomo
  stringMax=c(Kp_max, Ki_max, Kd_max),
  # Cromossomo inicial (nulo = aleatório)
  suggestions=NULL,
  # Tamanho da populaco
  popSize=populacao,
  # Quantidade de geracoes
  iters=geracoes,
  # Chance de mutacao ( padrao = 1/(size+1) )
  mutationChance=0.06,
  # Quantidade da populacao que passa para a proxima geracao
  elitism=NA,
  # Função de monitoramento de cada geracao
  monitorFunc=NULL,
  # Funcao de calculo de fitness
  evalFunc=fitness,
  # Mostra configuracoes
  showSettings=FALSE,
  # Escreve mais texto na tela
  verbose=FALSE
  )

# Obtemm os componentes do melhor cromossomo encontrado
bestKp <- GAmodel$population[which.min(GAmodel$evaluations),1]
bestKi <- GAmodel$population[which.min(GAmodel$evaluations),2]
bestKd <- GAmodel$population[which.min(GAmodel$evaluations),3]

# Imprime na tela o melhor cromossomo encontrado
cat("Melhor Kp:", bestKp, "\n")
cat("Melhor Ki:", bestKi, "\n")
cat("Melhor Kd:", bestKd, "\n")

# Mostra o grafico do melhor controlador obtido
pid = function(s) (bestKp + (bestKi/s) + bestKd*s)
malha_fechada = function(s) (pid(s)*sistema(s))/(1 + pid(s)*sistema(s))
resposta_degrau_unitario_f = function(s) (1/s)*malha_fechada(s)
resposta_do_sistema_t = invlap(Fs=resposta_degrau_unitario_f, 0, 2*pi, tx_amostragem)
plot(resposta_do_sistema_t, type="l")
erro = 0

for (i in 1:(tx_amostragem-1)){
  y_desejado = (1-exp(-10*resposta_do_sistema_t$x[i]))
  erro = erro + (y_desejado-resposta_do_sistema_t$y[i])^2
}
cat("Erro quadr�tico:", erro, "\n")

