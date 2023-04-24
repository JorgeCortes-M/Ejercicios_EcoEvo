#######################################################
########## Crecimiento Denso-dependiente ##############
#######################################################

#Para simular y graficar utilizaremos los siguientes paquetes

##install.packages("ggplot2")
##install.packages("ggpubr")
##instal.packages("latice")

# Nombre de los paquetes a utilizar
packages <- c("ggplot2", "zoo", "primer",
              "ggpubr", "lattice", "phaseR", "popbio")

# Instalar paquetes
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


library(ggplot2)
library(ggpubr)
library(primer)
library(zoo)
library(lattice)
library(phaseR)
library(popbio)

## podemos generar una función simple (en tiempo discreto) para 
## proyectar el tamaño poblacional futuro incorporando los parámetros del modelo logístico.

dlogistic <- function(K = 100, r = 1, N0 = 2, t = 15) {
    N <- c(N0, numeric(t))
      for (i in 1:t) N[i + 1] <- {
          N[i] + r * N[i] * (1 -  N[i]/K)
      }
          return(N)
}

## Podemos reemplazar los parámetros con valores predeterminados para 
## generar una proyección en el tiempo del tamaño poblacional y graficarlos: 

Nts <- dlogistic()

t <- 15; K <- 100
crec_logistic <- qplot(0:t, Nts, 
                       xlab = "Tiempo", 
                       ylab = "N",
                       main = "Crecimiento poblacional Denso-Dependiente", 
                       colour=I("Red"))
crec_logistic

## Incremento poblacional percápita vs Tamaño poblacional
## Exploremos como cambia el incremento total de la población en cada tiempo y el incremento percápita.

total.incr <- Nts[1:t + 1] - Nts[1:t]
per.capita.incr <- total.incr/Nts[1:t]


inc_total <- qplot(Nts[1:t], total.incr,
                   xlab = "N",
                   ylab = "Incremento total",
                   main = "Incremento total vs Tamaño poblacional",
                   colour=I("Red"))
inc_percapita <- qplot(Nts[1:t], per.capita.incr,
                       xlab = "N",
                       ylab = "Incremento percapita",
                       main = "Incremento percapita vs Tamaño poblacional",
                       colour=I("Red"))

inc_total
inc_percapita

ggarrange(crec_logistic, inc_percapita, inc_total)

## Análisis de estabilidad

## Podemos buscar los puntos de equilibrio utilizando el paquete phaseR

logisticPP <- phasePortrait(logistic, ylim = c(-5, 105), parameters = c(1, 100),
                            ylab="Tasa de crecimiento poblacional (dN/dt)", xlab="N")
N <- c(0, Nts[1:11])
pop.growth.rate <- expression( N * (1 -  N/100) )
points(N, eval(pop.growth.rate), cex=1.5)
text(N, eval(pop.growth.rate), letters[1:12],adj=c(.5,2))



## Efecto del tamaño inicial de la población

## Tomaremos 5 valores de población inicial y aplicaremos 
## nuestro modelo logístico a todos los valores para N0.

N0s <- c(0, 10, 20, 30, 150)
N <- data.frame(sapply(N0s, function(n) dlogistic(N0 = n)))

names(N)[names(N) == "X1"] <- "0"
names(N)[names(N) == "X2"] <- "10"
names(N)[names(N) == "X3"] <- "20"
names(N)[names(N) == "X4"] <- "30"
names(N)[names(N) == "X5"] <- "150"

autoplot(zoo(N), facet = NULL) + 
  geom_point() + 
  annotate(geom="text", x=13, y=20, label= "K=1/a", color = "red")+
  ggtitle(label = "Efecto del tamaño inicial de la población")+
  xlab(label = "Tiempo")+
  ylab(label = "N")


## ¿Qué efecto puedes notar al cambiar el tamaño inicial de la población?


## Efecto del valor de K

## Tomaremos 4 valores de K y le aplicaremos nuestro modelo logístico.

K.s <- c(150, 200, 600, 800)
N <- data.frame(sapply(K.s, function(a) dlogistic(K = a, t = 15)))

names(N)[names(N) == "X1"] <- "150"
names(N)[names(N) == "X2"] <- "200"
names(N)[names(N) == "X3"] <- "600"
names(N)[names(N) == "X4"] <- "800"

autoplot(zoo(N), facet = NULL) + geom_point() + 
  ggtitle(label = "Efecto de K sobre crecimiento Log")+
  xlab(label = "Tiempo")+
  ylab(label = "N")

## ¿Qué efecto puedes notar al cambiar los valores de K?


## Efecto de lambda sobre el crecimiento poblacional

## Tomaremos valores entre 1.3 y 2.8 para lambda de manera secuencial

lambda <- seq(1.3, 2.8, by = 0.3)
t <- 15
Ns <- data.frame(sapply(lambda, function(r) dlogistic(r = r,t = t)))

Ns

names(Ns)[names(Ns) == "X1"] <- "1.3"
names(Ns)[names(Ns) == "X2"] <- "1.6"
names(Ns)[names(Ns) == "X3"] <- "1.9"
names(Ns)[names(Ns) == "X4"] <- "2.2"
names(Ns)[names(Ns) == "X5"] <- "2.5"
names(Ns)[names(Ns) == "X6"] <- "2.8"

autoplot(zoo(Ns), facet = NULL) + geom_point() + 
  ggtitle(label = "Efecto de Lambda sobre crecimiento Log")+
  xlab(label = "Tiempo")+
  ylab(label = "N")

## ¿Qué efecto puedes notar al cambiar el valor de lambda? ¿Era algo esperado?

## Grafiquemos cada curva por separado para notar mejor el efecto de cambiar lambda.

tmp <- data.frame(lambda = as.factor(lambda), t(Ns))

Ns2 <- reshape(tmp, varying = list(2:ncol(tmp)), idvar = "lambda", v.names = "N", direction = "long")
str(Ns2)

print(xyplot(N ~ time | lambda, 
             data = Ns2, 
             type = "l", 
             layout = c(3, 2, 1), 
             col = 2,
             main="Curvas con distinto Lambda",
             xlab = "Tiempo"))

## ¿Qué puedes observar? ¿Existe algún patron?


## Anexo

## Diagramas de Bifurcación

## Ahora analizaremos qué ocurre al graficar N en función de distintos valores de lambda.

num.lambda <- 201; t <- 400
lambda.s <- seq(1, 3, length = num.lambda)

tmp <- sapply(lambda.s, function(r) dlogistic(r = r, N0 = 99, t = t))

tmp.s <- stack(as.data.frame(tmp))
names(tmp.s) <- c("N", "Old.Column.ID")
tmp.s$lambda <- rep(lambda.s, each = t + 1)
tmp.s$time <- rep(0:t, num.lambda)

## Tomamos sólo los valores finales que alcanza N para cada simulación.

N.bif <- subset(tmp.s, time > 0.5 * t)

d_bif <- qplot(lambda, N, data=N.bif,
      main = "Diagrama de Bifurcación",
      geom="point",
      size=I(0.5),
      colour=I("brown"))
d_bif+
  geom_vline(xintercept = 1.3, 
             linetype="dotted", 
            color = "red", 
            size=1.0) +
  geom_vline(xintercept = 1.6, 
             linetype="dotted", 
            color = "gold", 
            size=1.0) + 
  geom_vline(xintercept = 1.9, 
             linetype="dotted", 
             color = "green", 
             size=1.0)+ 
  geom_vline(xintercept = 2.2, 
             linetype="dotted", 
             color = "aquamarine", 
             size=1.0)+ 
  geom_vline(xintercept = 2.5, 
             linetype="dotted", 
             color = "blue", 
             size=1.0) + 
  geom_vline(xintercept = 2.8, 
             linetype="dotted", 
             color = "purple", 
             size=1.0)

## ¿Qué puedes observar? ¿Qué significan las bifurcaciones que se observan? ¿Cómo se relaciona este diagrama de bifurcación con lo observado en los gráficos anteriores?

## Diagramas de tela de araña

## Una forma de visualizar la estabilidad de un modelo poblacional es a través 
## de un diagrama de tela de araña, el cual considera el tamaño actual de la 
## población en el eje X y el tamaño en el tiempo siguiente en el eje Y.

out <- dlogistic(N0=2, r=2.9, K=100, t=1000)
n <- out[901:1000]
ntp1 <- n[-1]
nt <- n[-length(n)]
qplot(nt, ntp1, geom=c("point", "path") ) + labs(y="N[t+1]", x="N[t]")+geom_abline(intercept = 0, slope = 1, color = "red")

## ¿Qué observaríamos a tasas de crecimiento poblacional menores?


## Actividad

## Utiliza el set de datos Grizzly contenido en el paquete popbio, calcula los 
## valores de lambda, calcula un lambda promedio y utilizando un modelo 
## logístico con capacidad de carga 100 y como tamaño inicial el tamaño inicial 
## del set de datos Grizzly, genera los gráficos correspondientes del modelo 
## e indica si se genera una dinámica discreta normal, con oscilaciones o caótica. 
## Puedes repetir el proceso utilizando el valor de lambda más grande encontrado.

grizzly <- popbio::grizzly