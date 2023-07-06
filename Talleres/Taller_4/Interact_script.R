################# Interacciones Biológicas ############
#######################################################

# En este taller analizaremos 3 tipos de interacciones biológicas de forma
# resumida, nos centraremos en los modelos de Lotka-Volterra, pasando por simulaciones
# y en algunos casos intentaremos ajustar estos modelos a casos reales.

# Nombre de los paquetes a utilizar
packages <- c("tidyverse", "primer", "EcoVirtual", "FME", "phaseR", "devtools")

# Instalar paquetes
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Para iniciar cargaremos los paquetes a utilizar.

library(tidyverse)
library(EcoVirtual)
library(primer)
library(FME)
library(phaseR)
library(devtools)  

#Instalaremos además un paquete desde github

install_github("hallucigenia-sparsa/seqtime")  
library(seqtime) 



# La primera interacción a modelar será la competencia interespecífica
# Utilizaremos el modelo de Lotka Volterra

# Generaremos una función que aplica el modelo de Lotka-Volterra.
# Esta es una modificación de la que se puede encontrar en el paquete Ecovirtual


compLV=function(n01,n02,tmax,r1,r2,k1,k2,alfa,beta)
{
  resulta=matrix(0, ncol=3, nrow=tmax, dimnames=list(NULL, c("time", "Nsp1","Nsp2")))
  resulta[,1]=0:(tmax-1)
  resulta[1,c(2,3)]=c(n01,n02)
  for(t in 2:tmax)
  {
    nsp1=resulta[(t-1),2]
    nsp2=resulta[(t-1),3]
    resulta[t,2]=nsp1 + r1*nsp1*((k1-nsp1-alfa*nsp2)/k1)
    resulta[t,3]=nsp2 + r2*nsp2*((k2-nsp2-beta*nsp1)/k2)
    if (resulta[t,2]<1)  
    {
      resulta[t,2]=0
    }
    if (resulta[t,3]<1)  
    {
      resulta[t,3]=0
    }
  }
  old=par(mfrow=c(1,2), mar=c(4,4,2,1))
  plot(resulta[,1],resulta[,2],ylim=c(0,max(na.omit(resulta[,2:3]))),type="l",lty=4,xlab="time (t)",ylab="Population size", main="Population Growth", col="blue", lwd=1.5 )
  legend("topleft", legend=c("Sp. 1", "Sp. 2"), lty=4, col=c("blue", "green"), bty="n", cex=0.8)
  lines(resulta[,1],resulta[,3], col="green", lty=4, lwd=1.5)
  plot(resulta[,2],resulta[,3],type="l",col="red",xlab="N1",ylab="N2",ylim=c(0,max(c(na.omit(resulta[,3]),k1/alfa,k2))),xlim=c(0,max(c(na.omit(resulta[,2]),k2/beta,k1))), main="Isoclines")
  segments(0,k1/alfa,k1,0,lty=4, lwd=1.5, col="blue")
  segments(0,k2,k2/beta,0,lty=4,lwd=1.5, col="green" )
  legend("topleft", title="Equilibrium without habitat destruction",legend=c("isocline sp.1 ", "Isocline sp. 2", "Populations trajectory"), lty=c(4,4,1), col=c("blue", "green", "red"), bty="n", cex=0.8)
  invisible(resulta)
}

## Los parametros a ingresar son los tamaños inicales, las tasas de crecimiento,
## las capacidades de carga, los efectos de una especie sobre otra y el tiempo maximo.

## Partiremos de estos valores iniciales

compLV(n01=10, n02=10,r1=0.05, r2=0.03, k1=80, k2=50, alfa=1.2, beta=0.5, tmax=200)

## ¿A cual modelo visto en clases se parece?
## ¿Que ocurrira al aumentar levemente el tamaño inicial de una de las especies?

## Ahora probaremos el modelo con los siguientes parametros


compLV(n01=10, n02=10,r1=0.05, r2=0.03, k1=80, k2=50, alfa=2.0, beta=0.5, tmax=200)

## ¿A cual modelo visto en clases se parece?
## ¿Que especie es excluida competitivamente?

compLV(n01=10, n02=10,r1=0.05, r2=0.03, k1=80, k2=50, alfa=1.2, beta=0.8, tmax=200)

## ¿A cual modelo visto en clases se parece?
## ¿Que especie es excluida competitivamente?

compLV(n01=10, n02=10,r1=0.05, r2=0.03, k1=80, k2=50, alfa=3.0, beta=1.8, tmax=200)

## ¿A cual modelo visto en clases se parece?
## ¿Que especie es excluida competitivamente?
## Si modifican los tamaños iniciales ¿Cambia el resultado de la competencia?

#Ahora a modo de demostración veremos el modelo de predador presa


# parametros
pars <- c(alpha = 1, beta = 0.2, delta = 0.5, gamma = 0.2)
# estado inicial
init <- c(x = 1, y = 2)
# tiempo
times <- seq(0, 100, by = 1)

#Creamos la función que muestra las derivadas
deriv <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    d_x <- alpha * x - beta * x * y
    d_y <- delta * beta * x * y - gamma * y
    return(list(c(x = d_x, y = d_y)))
  })
}

# Resolvemos la ecuación diferencial usando los parametros iniciales y el tiempo

lv_results <- ode(init, times, deriv, pars)


# Graficamos

lv_results %>% 
  data.frame() %>% 
  gather(var, pop, -time) %>% 
  mutate(var = if_else(var == "x", "Presa", "Predador")) %>% 
  ggplot(aes(x = time, y = pop)) +
  geom_line(aes(color = var)) +
  scale_color_brewer(NULL, palette = "Set1") +
  labs(title = "Modelo de predador presa de Lotka-Volterra",
       subtitle = paste(names(pars), pars, sep = " = ", collapse = "; "),
       x = "Tiempo", y = "N")

# Generamos el diagrama de fase para estas ecuaciones

lotkaVolterra_flowField    <- flowField(lotkaVolterra,
                                        xlim       = c(0, 10), 
                                        ylim       = c(0, 10),
                                        parameters = c(1, 0.2, 0.1, 0.2),
                                        add        = F,
                                        xlab = "N Presa",
                                        ylab = "N Predador")


lotkaVolterra_trajectories <- trajectory(lotkaVolterra,
                                         y0     = rbind(c(2, 2),
                                                        c(0.5, 0.5),
                                                        c(0.5, 1.5),
                                                        c(1.5, 0.5),
                                                        c(3, 3)),
                                         parameters = c(1, 0.2, 0.1, 0.2),
                                         col    = rep("black", 5),
                                         tlim   = c(0, 100))

lotkaVolterra_isoclines <- nullclines(lotkaVolterra,
                                         xlim      = c(0,10),
                                         ylim      = c(0,10),
                                         col = c("aquamarine2", "blueviolet"),
                                         parameters = c(1, 0.2, 0.1, 0.2),
                                         points = 251)


# Mutualismo, en este caso los efectos percapita de una especie sobre la otra
# son positivos


mutLV=function(n01,n02,tmax,r1,r2,k1,k2,alfa,beta)
{
  resulta=matrix(0, ncol=3, nrow=tmax, dimnames=list(NULL, c("time", "Nsp1","Nsp2")))
  resulta[,1]=0:(tmax-1)
  resulta[1,c(2,3)]=c(n01,n02)
  for(t in 2:tmax)
  {
    nsp1=resulta[(t-1),2]
    nsp2=resulta[(t-1),3]
    resulta[t,2]=nsp1 + nsp1*r1*((k1 - nsp1 + alfa*nsp2))/k1
    resulta[t,3]=nsp2 + nsp2*r2*((k2 - nsp2 + beta*nsp1))/k2
    if (resulta[t,2]<1)  
    {
      resulta[t,2]=0
    }
    if (resulta[t,3]<1)  
    {
      resulta[t,3]=0
    }
  }
 
  old=par(mfrow=c(1,2), mar=c(4,4,2,1))
  plot(resulta[,1],resulta[,2],ylim=c(0,max(na.omit(resulta[,2:3]))),type="l",lty=4,xlab="time (t)",ylab="Population size", main="Population Growth", col="blue", lwd=1.5 )
  legend("topleft", legend=c("Sp. 1", "Sp. 2"), lty=4, col=c("blue", "green"), bty="n", cex=0.8)
  lines(resulta[,1],resulta[,3], col="green", lty=4, lwd=1.5)
  plot(resulta[,2],resulta[,3],type="l",col="red",xlab="N1",ylab="N2",ylim=c(0,max(c(na.omit(resulta[,3]),-k1/alfa,k2))),xlim=c(0,max(c(na.omit(resulta[,2]),-k2/beta,k1))), main="Isoclines")
  segments(k1,0,k1+alfa*max(resulta[,3]),max(resulta[,3]),lty=4, lwd=1.5, col="blue")
  segments(0,k2,max(resulta[,2]),k2+max(resulta[,2])*beta,lty=4,lwd=1.5, col="green" )
  legend("topleft", title="Equilibrium without habitat destruction",legend=c("isocline sp.1 ", "Isocline sp. 2", "Populations trajectory"), lty=c(4,4,1), col=c("blue", "green", "red"), bty="n", cex=0.8)
  invisible(resulta)
}


mutLV(n01 = 10,n02 = 10,tmax = 300,r1 = 0.01,r2 = 0.02,k1 = 100,k2 = 150,alfa = 0.9,beta = 0.8)

# Cambiando los parámetros intenta generar las dinámicas vistas en clase

## Ejemplo de lotka volterra generalizado

S <- 10
alpha <- .01
r <- runif(S)*2
a <- matrix(rnorm(S^2, m=alpha, sd=alpha/10), nrow=S, ncol=S)
parms <- list(r,a)
t=seq(0,40, by=.1)
N0 <- runif(S)/(S*alpha)
library(deSolve)
lvout <- ode(N0, t, lvcompg, parms)
matplot(t, lvout[,-1], type="l", ylab="N", log='y')


# Una forma más simple

tsplot(glv(N=10,generateA(10)),header="gLV")
