####### Configurações do programa #######

# Modelo do sistema de segunda ordem (as variáveis importantes)
alpha = 1
wn = 100
csi = 0.5

# Configurações do algoritmo genético
geracoes = 100
populacao = 10

Kp_min = 0
Kp_max = 30

Ki_min = 0
Ki_max = 100

Kd_min = 0
Kd_max = 100

#########################################

# Biblioteca que contém invlap() - transformada inversa de laplace
library(pracma)
library(genalg)

# Função que avalia o cromossomo
# cromossomo = (P, I, D)
fitness = function(cromossomo) {
  Kp = cromossomo[1]
  Ki = cromossomo[2]
  Kd = cromossomo[3]
  
  # Esta função de teste busca otimizar Kp para que seja igual à 20.
  # Ki e Kd não estão sendo avaliados.
  return(abs(Kp-20))
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
cat("Melhor Kp: ", bestKp, "\n")
cat("Melhor Ki: ", bestKi, "\n")
cat("Melhor Kd: ", bestKd, "\n")