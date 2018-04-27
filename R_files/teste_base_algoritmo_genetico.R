####### Configurações do programa #######

# Modelo do sistema de segunda ordem (as variáveis importantes)
alpha = 1
wn = 1000
csi = 0.5

# Configurações do algoritmo genético
geracoes = 100
populacao = 10

Kp_min = 0
Kp_max = 50

Ki_min = 0
Ki_max = 50

Kd_min = 0
Kd_max = 50

#########################################

# Modelo do sistema de segunda ordem
sistema = function(s) (alpha*(wn^2))/(s^2 + 2*csi*wn*s + wn^2)

# Biblioteca que contém invlap() - transformada inversa de laplace
library(pracma)
library(genalg)

# Função que avalia o cromossomo
# cromossomo = (P, I, D)
fitness = function(cromossomo) {
  # Separa o cromossomo nas variáveis de ajuste do PID
  Kp = cromossomo[1]
  Ki = cromossomo[2]
  Kd = cromossomo[3]
  
  # Modelo do controlador PID
  pid = function(s) (Kp + (Ki/s) + Kd*s)
  
  # PID e sistema com realimentação recebendo entrada de um degrau
  malha_fechada = function(s) (pid(s)*sistema(s))/(1 + pid(s)*sistema(s))

  # Resposta do sistema realimentado ao degrau unitário no domínio da frequência
  resposta_degrau_unitario_f = function(s) (1/s)*malha_fechada(s)
  
  # Resposta do sistema com controlador ao degrau unitário no domínio do tempo
  resposta_do_sistema_t = invlap(Fs=resposta_degrau_unitario_f, t1=0, t2=2*pi, nnt=1000)
  
  # A função de fitness retorna a soma do erro quadrático em relação à um degrau unitário (resposta utópica)
  erro = 0
  for (y in resposta_do_sistema_t$y){
    erro = erro + (1-y)^2
  }
  return(erro)
}

# Algoritmo genético (https://www.rdocumentation.org/packages/genalg/versions/0.2.0/topics/rbga)
GAmodel = rbga(
  # Valores mínimos de cada componente do cromossomo
  stringMin=c(Kp_min, Ki_min, Kd_min),
  # Valores máximos de cada componente do cromossomo
  stringMax=c(Kp_max, Ki_max, Kd_max),
  # Cromossomo inicial (nulo = aleatório)
  suggestions=NULL,
  # Tamanho da população
  popSize=populacao,
  # Quantidade de gerações
  iters=geracoes,
  # Chance de mutação ( padrão = 1/(size+1) )
  mutationChance=NA,
  # Quantidade da população que passa para a próxima geração
  elitism=NA,
  # Função de monitoramento de cada geração
  monitorFunc=NULL,
  # Função de cálculo de fitness
  evalFunc=fitness,
  # Mostra configurações
  showSettings=FALSE,
  # Escreve mais texto na tela
  verbose=FALSE
  )

# Obtém os componentes do melhor cromossomo encontrado
bestKp <- GAmodel$population[which.min(GAmodel$evaluations),1]
bestKi <- GAmodel$population[which.min(GAmodel$evaluations),2]
bestKd <- GAmodel$population[which.min(GAmodel$evaluations),3]

# Imprime na tela o melhor cromossomo encontrado
cat("Melhor Kp:", bestKp, "\n")
cat("Melhor Ki:", bestKi, "\n")
cat("Melhor Kd:", bestKd, "\n")

# Mostra o gráfico do melhor controlador obtido
pid = function(s) (bestKp + (bestKi/s) + bestKd*s)
malha_fechada = function(s) (pid(s)*sistema(s))/(1 + pid(s)*sistema(s))
resposta_degrau_unitario_f = function(s) (1/s)*malha_fechada(s)
resposta_do_sistema_t = invlap(Fs=resposta_degrau_unitario_f, t1=0, t2=4*pi, nnt=1000)
plot(resposta_do_sistema_t, type="l")
erro = 0
for (y in resposta_do_sistema_t$y){
  erro = erro + (1-y)^2
}
cat("Erro quadrático:", erro, "\n")

