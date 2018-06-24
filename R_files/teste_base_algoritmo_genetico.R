# Importa o csv com a resposta ao degrau

RespostaDegrau <- read.csv(file="resposta2.csv", header=TRUE, sep=",");

# Determina o overshoot da resposta

Mp = max(RespostaDegrau[1]);
indexPico = match(Mp, unlist(RespostaDegrau[1]));
Tp = RespostaDegrau[indexPico,2];
Mp = Mp-1

cat("Tp =",Tp,"\n")
cat("Mp =", Mp,"\n")


# Cálculo para determinação do csi
csi = sqrt((log(Mp)/pi)^2/(1+(log(Mp)/pi)^2))

cat("csi =",csi,"\n")

# Cálculo para determinação do wn
wn = pi/(Tp*sqrt(1-csi^2))

cat("wn =", wn,"\n")

# Determinacao do periodo de amostragem
Ta = Tp

# Determinacao da funcao de comparacao
Tc = Tp*100

# Configuracoes do algoritmo genetico
geracoes = 50
populacao = 100
mutacao = 0.06


Kp_min = 0
Kp_max = 100

Ki_min = 0
Ki_max = 100

Kd_min = 0
Kd_max = 100

tx_amostragem = 5000

#valores iniciais
init = c(1,0.1,0.1)
init_df = data.frame(init)

#########################################

# Modelo do sistema de segunda ordem
sistema = function(s) ((wn^2))/(s^2 + 2*csi*wn*s + wn^2)

# Biblioteca que contem invlap() - transformada inversa de laplace
library(pracma)
library(genalg)

erro = 0
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
  
  # Resposta do sistema com controlador ao degrau unitÃ¡rio no dominio do tempo
  resposta_do_sistema_t = invlap(resposta_degrau_unitario_f, 0, Ta, tx_amostragem)
  
  # A funcao de fitness retorna a soma do erro quadratico em relacao a uma exponencial (resposta utopica)
  erro = 0
  for (i in 1:(tx_amostragem-1)){
    y_desejado = (1-exp(-Tc*resposta_do_sistema_t$x[i]))
    # y_desejado = 1
    erro = erro + (y_desejado-resposta_do_sistema_t$y[i])^2
  }
  return(erro)
}

# funcao de monitoramenteo, indica os melhores valores para pid em cada geracao indicada
monitor <- function(obj) {

if (obj$iter == 1) {
  print(paste("GENERATION:", obj$iter))
  print(obj$population[which.max(obj$evaluations), ])
}
if (obj$iter == 10) {
  print(paste("GENERATION:", obj$iter))
  print(obj$population[which.max(obj$evaluations), ])
}
if (obj$iter == 20) {
  print(paste("GENERATION:", obj$iter))
  print(obj$population[which.min(obj$evaluations), ])
}
if (obj$iter == 30) {
  print(paste("GENERATION:", obj$iter))
  print(obj$population[which.min(obj$evaluations), ])
}
if (obj$iter == 40) {
  print(paste("GENERATION:", obj$iter))
  print(obj$population[which.min(obj$evaluations), ])
}

}

# Algoritmo genetico (https://www.rdocumentation.org/packages/genalg/versions/0.2.0/topics/rbga)
GAmodel = rbga(
  # Valores máxnimos de cada componente do cromossomo
  stringMin=c(Kp_min, Ki_min, Kd_min),
  # Valores maximos de cada componente do cromossomo
  stringMax=c(Kp_max, Ki_max, Kd_max),
  # Cromossomo inicial (nulo = aleatorio)
  suggestions= init_df,
  # Tamanho da populaco
  popSize=populacao,
  # Quantidade de geracoes
  iters=geracoes,
  # Chance de mutacao ( padrao = 1/(size+1) )
  mutationChance=mutacao,
  # Quantidade da populacao que passa para a proxima geracao
  elitism=NA,
  # FunÃ§Ã£o de monitoramento de cada geracao
  monitorFunc=monitor,
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
resposta_do_sistema_t = invlap(Fs=resposta_degrau_unitario_f, 0, Ta, tx_amostragem)

x = resposta_do_sistema_t$x
matplot(x, cbind(resposta_do_sistema_t$y,1-exp(-Tc*x)), type = "l", col = c("blue","red"), ylim = c(0,1.3))
      

erro = 0
for (i in 1:(tx_amostragem-1)){
  y_desejado = (1-exp(-Tc*resposta_do_sistema_t$x[i]))
  # y_desejado = 1
  erro = erro + (y_desejado-resposta_do_sistema_t$y[i])^2
}

cat("Erro quadrático:", erro, "\n")

