################################################
### CRECIMIENTO DENSO-INDEPENDIENTE DISCRETO ###
################################################

##Crecimiento poblacional de Nymphaea odorata a través de los años##

#install.packages("primer")
#install.packages("ggplot2")
#install.packages("ggpubr")

# Nombre de los paquetes a utilizar
packages <- c("ggplot2", "primer",
              "ggpubr", "EcoVirtual", "popbio")

# Instalar paquetes
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


library(primer)
library(ggplot2)
library(ggpubr)
library(popbio)


N <- c(1, 3, 9, 27, 81)
year <- 2001:2005
qplot(year, N, 
      xlab = "Año", 
      ylab = "N", 
      main = "Conteo por año de N. odorata", 
      colour=I("Red"))
##calculamos las tasas de crecimiento dividiendo cada número deindividuos de un año por el número
##de inviduos del año anterior: lambda = Nt+1/Nt

rates = N[2:5]/N[1:4] #divide cada valor de N por el valor anterior
rates

## teniendo esto en cuenta podríamos escribir un modelo denso independiente sencillo
## Nt = lambda^t * N0, para este caso lambda es igual a 3. (tasa finita de crecimiento)
## Damos parametros iniciales a nuestro modelo

N0 <- 1
lambda <- 2
time <- 0:10

Nt <- N0*lambda^time

qplot(time, Nt, 
      main = "Modelo denso independiente discreto",
      xlab = "Tiempo", 
      ylab = "N", 
      colour=I("Blue"))

## Exploremos el efecto del tamaño inicial de la población

N0  <- c(10,20,30)
lambda <- 2
time <- 0:4

##utilizaremos sapply para aplicar nuestro modelo a cada elemento en N0, es decir a cada
##tamaño poblacional inicial

Nt.s <- sapply(N0, function(n) n * lambda^time)
Nt.s

##usamos la función matplot para graficar las 3 curvas a la vez
##aplicamos logaritmo a los tama?os poblacionales.

par(mfrow=c(1,2))

matplot(time, Nt.s, type = "l", 
        pch = 1:3, 
        xlab = "Tiempo", 
        ylab = "N", 
        main="Tamaño poblacional vs Tiempo")

matplot(time, Nt.s, type = "l",
        log = "y", 
        pch = 1:3, 
        xlab = "Tiempo", 
        ylab = "log(N)", 
        main="Log Tamaño poblacional vs Tiempo")

##¿Qué efectos puedes notar al cambiar la población inicial?

##Efecto del valor de lambda

## Utilizamos distintos valores de lambda y vemos su efecto sobre el crecimiento poblacional

N0 <- 100
time <- 0:3
lambdas <- c(0.5, 1, 1.5)


##nuevamente utilizamos sapply para aplicar nuestro modelo denso independiente a todos los
## valores de lambda

N.all <- sapply(lambdas, function(x) N0 * x^time)

par(mfrow=c(1,2))

matplot(time, N.all, type = "l", 
        main = "Crecimiento exp. continuo",
        xlab = "Tiempo",
        ylab = "N", 
        col = 6)
legend("topleft", paste(rev(lambdas)), lty = 5:1, col = 1, bty = "n", title = "lambdas")
matplot(time, N.all, type = "l", 
        main = "Log(N) vs Tiempo",
        xlab = "Tiempo", 
        ylab = "log(N)", 
        log = "y", 
        col = 6)


## BOX1
## Simulaciones de Crecimiento exponencial en el paquete EcoVirtual


N0 = 10
lamb= 1.5
tmax= 10
par(mfrow=c(1,1))
popExp(N0, lamb, tmax, intt = 1)


## Otra métrica que podemos explorar además de r para analizar el crecimiento de una población
## es el tiempo en que se duplica su población, el cual se obtiene de:
##                    t = ln(2)/r

## Creamos una función para explorar esto:

m.time <- function(r, m = 2) {
  log(m)/r
}

rs <- c(0, 0.1, 0.5, 1, 2)
times <- m.time(rs)

qplot(x = rs, y = times,
      colour = "red",
      size = 2.0)+
  geom_line(size=1.0)+
  ggtitle("Tiempo en que se duplica una población en función de r")+
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    )

##ANEXO

################################################
### CRECIMIENTO DENSO-INDEPENDIENTE CONTINUO ###
################################################

##Actividad

## El crecimiento denso-independiente o exponenial continuo puede ser modelado por la ecuaci?n:
##                       Nt = N0e^(rt)
## Esta ecuaci?n puede ser obtenida de la resolución de un límite (ver en Anexo)
## Podemos proyectar en el tiempo la población tal como hicimos con el modelo en su forma discreta:

## Escogemos 5 tasas de crecimiento (r) y un intervalo de tiempo de 1 a 100

r <- c(-0.03, -0.02, 0, 0.02, 0.03)
N0 <- 2; t <- 1:100
cont.mat <- sapply(r, function(ri) N0 * exp(ri * t))

## Graficamos

layout(matrix(1:2, nrow = 1))
matplot(t, cont.mat, type = "l", 
        main = "Crecimiento exp. continuo",
        xlab = "Tiempo",
        ylab = "N", 
        col = 6)
legend("topleft", paste(rev(r)), lty = 5:1, col = 1, bty = "n", title = "r")

matplot(t, cont.mat, type = "l", 
        main = "Log(N) vs Tiempo",
        xlab = "Tiempo", 
        ylab = "log(N)", 
        log = "y", 
        col = 6)



## Hasta ahora hemos utilizado tasas fijas de lambda, pero este valor podría cambiar en cada año
## podríamos tomar por ejemplo la media aritmética de los distintos valores de lambda
## o podríamos tomar la media geométrica, a continuación las comparamos.
## usaremos el set de datos sparrows contenidos en el paquete primer

data(sparrows)
str(sparrows)
attach(sparrows)
qplot(x = Year, y = Count, data = sparrows, 
      xlab = "A?o",
      ylab = "Conteo",
      geom = "point", 
      colour=I("Blue"),
      main = "Conteo de Melospiza melodia")


##usaremos 6 valores observados de R de los datos sparrows

t <- 5
SS6 <- sparrows[1:(t + 1), ]

#Calculamos los valores de lambda para cada intervalo y calculamos media aritm?tica y geom?trica.

SSgr <- SS6$Count[2:(t + 1)]/SS6$Count[1:t]
lam.A <- sum(SSgr)/t
lam.G <- prod(SSgr)^(1/t)

lam.A
lam.G

#Ahora graficamos las proyecciones de ambas medias

N0 <- SS6$Count[1]
Nt_lamA <- (sapply(N0, function(n) n * lam.A^(0:t)))
Nt_lamG <- (sapply(N0, function(n) n * lam.G^(0:t)))

data_fin <- cbind(SS6, Nt_lamA, Nt_lamG)

colors <- c("Nt_lamA" = "red", "Nt_lamG" = "blue")

ggplot(data_fin, aes(x = Year, y = Count)) +
  ggtitle("Población proyectada")+
  geom_point(size = 4, shape = 20, color = "orange")+
  geom_line(data = data_fin, aes(x=Year, y=Nt_lamA, color = "Nt_lamA"), linetype = "dashed")+
  geom_line(data = data_fin, aes(x=Year, y=Nt_lamG, color = "Nt_lamG"))+
  labs(x = "Años",
       y = "N (individuos)",
       color="Legend")+
  scale_color_manual(values = colors)

##A partir de estos resultados, averigua en qué se diferencia la media aritmética y la media geométrica.

## Actividad
## A continuación se encuentran los datos de crecimiento poblacional de
## Osos Grizzly en Yellowstone, grafica la dinámica poblacional, identifica
## si hay un segmento donde pudo haber crecimiento exponencial
## Calcula los valores de lambda de dicho segmento
## Calcula un valor promedio de lambda y simula a 20 años el crecimiento poblacional

grizzly <- popbio::grizzly

