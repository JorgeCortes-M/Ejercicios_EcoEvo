## Estocasticidad Ambiental ##


## Cargamos los paquetes necesarios

## install.packages("ggplot2")
## install.packages("primer")
## install.packages("EcoVirtual") 
## install.packages("fitdistrplus")
## install.packages("EnvStats")
## install.packages("remotes")

remotes::install_github("BruceKendall/PVA")


# Nombre de los paquetes a utilizar
packages <- c("ggplot2", "primer", "EcoVirtual", "EnvStats", "fitdistrplus")

# Instalar paquetes
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


library(ggplot2)
library(primer)
library(EcoVirtual)
library(EnvStats)
library(fitdistrplus)
library(PVA)

# Cargamos los datos de Melospiza melodia contenidos en el paquete primer

data(sparrows)
names(sparrows)
attach(sparrows)

## Graficamos los conteos en el tiempo de la poblaci?n.

qplot(Year, Count,
      geom = c("point", "line"),
      xlab = "Años",
      ylab = "Conteo",
      main = "Conteo de M. melodia",
      colour=I("red"))+ theme_bw()


## Calculamos lambda y graficamos.

obs.lam <- Count[-1]/Count[-length(Count)]

qplot(Year[-length(Count)], obs.lam,
      geom = c("point", "line"),
      xlab = "Años",
      ylab = "Tasa de crecimiento poblacional",
      main = "Tasa de crecimiento de M. melodia",
      colour = I("brown"))+
  geom_hline(yintercept=1, linetype="dashed", color="red", size = 1)+
  theme_bw()

## Graficamos el histograma de frecuencias de lambda.




hist(obs.lam, 
     xlim = c(0,5),
     ylim = c(0,10), 
     breaks = 20,
     ylab = "Frecuencia",
     xlab = "Lambda",
     main = "Histograma de frecuencia de lambda",
     col = "orange")

## Intentamos identificar una distribución de estos datos

descdist(obs.lam, discrete = FALSE)
fit.lognorm <- fitdist(obs.lam, "lnorm")
fit.gamma <- fitdist(obs.lam, "gamma")

plot(fit.lognorm)
plot(fit.gamma)

fit.lognorm$aic
fit.gamma$aic


## Calculamos algunas estadísticas de lambda observado



lam_stats <- elnormAlt(obs.lam, method = "mvue", ci = T, ci.type = "two-sided",
          ci.method = "land", conf.level = 0.95, parkin.list = NULL)

lam_stats$parameters

lam_stats$interval$limits

##PVA##

# Hacemos una realización para la población a 100 años y graficamos.

years <- 100
set.seed(3)
sim.Rs <- sample(x = obs.lam, size = years, replace = TRUE)
output <- numeric(years + 1)

output[1] <- Count[Year == max(Year)]

for (t in 1:years) output[t + 1] <- {
  
  output[t] * sim.Rs[t]
}

par(mfrow=c(1,1))
Ye <- 0:100
qplot(Ye, output, geom = c("point","line")
      ,colour = I("red"),
      xlab = "Años",
      ylab = "Tamaño poblacional",
      main = "1 realización")+
  theme_bw()

## Haremos 10 realizaciones a partir de la población original y graficamos.

sims = 10
sim.RM <- matrix(sample(obs.lam, sims * years, replace = TRUE),
                 
                 nrow = years, ncol = sims)


output[1] <- Count[Year == max(Year)]
outmat <- sapply(1:sims, function(i) {
  
  for (t in 1:years) output[t + 1] <- output[t] * sim.RM[t,i]
  
  output
})

par(mfrow=c(1,1))

matplot(0:years, outmat, type = "l", 
        log = "y", 
        main = "10 realizaciones",
        ylab = "log(N)",
        xlab = "Años")



## Generamos una funci?n para generar todas las realizaciones que deseemos

PopSim <- function(Rs, N0, years = 100, sims = 10) {
  
  sim.RM = matrix(sample(Rs, size = sims * years, replace = TRUE),
                  
                  nrow = years, ncol = sims)
  
  output <- numeric(years + 1)
  
  output[1] <- N0
  
  outmat <- sapply(1:sims, function(i) {
    
    for (t in 1:years) output[t + 1] <- round(output[t] *
                                                
                                                sim.RM[t, i], 0)
    
    output
    
  })
  
  return(outmat)
}

## Probamos la funci?n con los valores de lambda calculados y 1000 realizaciones

output <- PopSim(Rs = obs.lam, N0 = 43, sims = 1000, years = 500)

N.500 <- output[501, ]
summary(N.500, digits = 6)

## Graficamos los tamaños poblacionales a 100 años de todas las realizaciones
N.100 <- output[101, ]
hist(N.100, breaks = 20000,
     xlim=c(1,10000))

## Graficamos las realizaciones generadas

matplot(0:500, output, type = "l", 
        log = "y", 
        main = "1000 realizaciones",
        ylab = "log(N)",
        xlab = "Años")


## Generamos una función para obtener las probabilidades de extinción en cada año y graficamos.

Ext <- function(output, sims=1000){
  num_ex <- 0
  prob_ex <- 0
  for(i in 0:501){
    num_ex[i] <- sum(ifelse(output[i,] == 0, 1, 0))
    prob_ex[i] <- num_ex[i]/sims
  }
  return(prob_ex)
}

proba <- Ext(output,sims=1000)

qplot(0:500, proba,
      geom = "line",
      xlab = "Años",
      ylab = "Probabilidad de extinción",
      main = "Curva de quasi-extinción de M. melodia",
      colour = I("red"))+
  theme_bw()

## Paquete PVA

melodia.ext <- count.DI.PVA(Count, Nx=1, nboot=1000, year.max=500)

melodia.ext


## Ahora es tu turno de realizar este análisis con los datos de Grizzly que se encuentran a continuaci?n.
## Genera las estadísticas y gráficos que caracterizan a los datos y un análisis de PVA.

## Secuencia de años

year=seq(from=1959,to=1996,by=1)
##########

# Conteo

N=c(44,47,46,44,46,45,46,40,39,39,42,39,41,40,33,36,34,39,35,34,38,36,37,41,39,51,
    47,57,48,60,65,74,69,65,57,70,81,99)


lam <-  N[-1]/N[-length(N)]

output <- PopSim(Rs = lam, N0 = 10, sims = 1000, years = 500)

proba <- Ext(output,sims=1000)

qplot(year, N,
      geom = "line",
      xlab = "Años",
      ylab = "Probabilidad de extinción",
      main = "Curva de quasi-extinción de M. melodia",
      colour = I("red"))+
  theme_bw()

qplot(0:500, proba,
      geom = "line",
      xlab = "Años",
      ylab = "Probabilidad de extinción",
      main = "Curva de quasi-extinción de M. melodia",
      colour = I("red"))+
  theme_bw()

melodia.ext <- count.DI.PVA(N, Nx=60, nboot=1000, year.max=500)

qplot(year[-length(Count)], lam,
      geom = c("point", "line"),
      xlab = "Años",
      ylab = "Tasa de crecimiento poblacional",
      main = "Tasa de crecimiento de M. melodia",
      colour = I("brown"))+
  geom_hline(yintercept=1, linetype="dashed", color="red", size = 1)+
  theme_bw()

PVA::