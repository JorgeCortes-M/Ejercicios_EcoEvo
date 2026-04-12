#######################################################
########## Crecimiento Denso-dependiente ##############
#######################################################

# Paquetes necesarios para el taller
packages <- c("ggplot2", "ggpubr", "lattice", "phaseR", "popbio", "deSolve")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

library(ggplot2)
library(ggpubr)
library(lattice)
library(phaseR)
library(popbio)
library(deSolve)

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

ggplot(data.frame(time = 0:t, N = Nts), aes(x = time, y = N)) +
  geom_line(color = "red", linewidth = 1.2) +
  geom_point(color = "red", size = 2) +
  geom_hline(yintercept = K, linetype = "dashed", color = "gray50") +
  labs(x = "Tiempo", y = "N", title = "Crecimiento poblacional Denso-Dependiente") +
  theme_bw()

## Incremento poblacional percápita vs Tamaño poblacional
## Exploremos como cambia el incremento total de la población en cada tiempo y el incremento percápita.

total.incr <- Nts[1:t + 1] - Nts[1:t]
per.capita.incr <- total.incr/Nts[1:t]

df_incr <- data.frame(
  N = Nts[1:t],
  total = total.incr,
  per_capita = per.capita.incr
)

p1 <- ggplot(df_incr, aes(x = N, y = total)) +
  geom_point(color = "red", size = 2) +
  geom_line(color = "red") +
  labs(x = "N", y = "Incremento total", title = "Incremento total vs Tamaño poblacional") +
  theme_bw()

p2 <- ggplot(df_incr, aes(x = N, y = per_capita)) +
  geom_point(color = "blue", size = 2) +
  geom_line(color = "blue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(x = "N", y = "Incremento per cápita", title = "Incremento per cápita vs Tamaño poblacional") +
  theme_bw()

p3 <- ggplot(data.frame(time = 0:t, N = Nts), aes(x = time, y = N)) +
  geom_line(color = "darkgreen", linewidth = 1.2) +
  geom_point(color = "darkgreen", size = 2) +
  geom_hline(yintercept = K, linetype = "dashed", color = "gray50") +
  labs(x = "Tiempo", y = "N", title = "Crecimiento logístico") +
  theme_bw()

ggarrange(p3, p2, p1, ncol = 3)

## Análisis de estabilidad

## Los puntos de equilibrio son N = 0 y N = K
## Visualicemos el diagrama de fase

N_seq <- seq(0, 150, length = 100)
r <- 1
K <- 100
dNdt <- r * N_seq * (1 - N_seq / K)

ggplot(data.frame(N = N_seq, dNdt = dNdt), aes(x = N, y = dNdt)) +
  geom_line(color = "purple", linewidth = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = K, linetype = "dashed", color = "red") +
  annotate("text", x = K + 5, y = 25, label = "K = 100", color = "red") +
  labs(x = "N (tamaño poblacional)", y = "dN/dt", title = "Diagrama de fase: Modelo logístico") +
  theme_bw()

## Efecto del tamaño inicial de la población

## Tomaremos 5 valores de población inicial y aplicaremos 
## nuestro modelo logístico a todos los valores para N0.

N0s <- c(0, 10, 20, 30, 150)
t <- 15

N_all <- sapply(N0s, function(n) dlogistic(N0 = n, t = t))

df_N0 <- data.frame(
  time = rep(0:t, length(N0s)),
  N = as.vector(N_all),
  N0 = rep(as.character(N0s), each = t + 1)
)

ggplot(df_N0, aes(x = time, y = N, color = N0, group = N0)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray50") +
  labs(x = "Tiempo", y = "N", title = "Efecto del tamaño inicial de la población", color = "N0") +
  theme_bw()

## Efecto del valor de K

## Tomaremos 4 valores de K y le aplicaremos nuestro modelo logístico.

K.s <- c(150, 200, 600, 800)
t <- 15

N_K <- sapply(K.s, function(k) dlogistic(K = k, t = t))

df_K <- data.frame(
  time = rep(0:t, length(K.s)),
  N = as.vector(N_K),
  K = rep(as.character(K.s), each = t + 1)
)

ggplot(df_K, aes(x = time, y = N, color = K, group = K)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  labs(x = "Tiempo", y = "N", title = "Efecto de K sobre el crecimiento logístico") +
  theme_bw()

## Pregunta: ¿Cómo cambia la forma de la curva al variar K?


## Efecto de la tasa de crecimiento (r)

## En el modelo discreto, la relación entre lambda y r es: lambda = 1 + r
## Veamos qué ocurre al variar la tasa de crecimiento:

r_values <- seq(0.1, 1.5, by = 0.2)
t <- 20

N_r <- sapply(r_values, function(r) dlogistic(r = r, K = 100, t = t))

df_r <- data.frame(
  time = rep(0:t, length(r_values)),
  N = as.vector(N_r),
  r = rep(round(r_values, 1), each = t + 1)
)

ggplot(df_r, aes(x = time, y = N, color = factor(r), group = factor(r))) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray50") +
  labs(x = "Tiempo", y = "N", title = "Efecto de r sobre el crecimiento logístico", color = "r") +
  theme_bw()

## Pregunta: Para valores grandes de r (> 1), ¿qué sucede con la dinámica poblacional?

################################################
### CRECIMIENTO LOGÍSTICO CONTINUO ###
################################################

## El modelo logístico continuo (ecuación diferencial) es:
## dN/dt = rN * (1 - N/K)
## Donde:
##   r = tasa de crecimiento intrínseca
##   K = capacidad de carga

## La solución analítica de esta EDO es:
## N(t) = K / (1 + ((K - N0) / N0) * e^(-r*t))

logistic_solution <- function(t, N0, r, K) {
  K / (1 + ((K - N0) / N0) * exp(-r * t))
}

t <- 0:50
N0 <- 10
r <- 0.5
K <- 100

N_logistic <- logistic_solution(t, N0, r, K)

ggplot(data.frame(t, N_logistic), aes(x = t, y = N_logistic)) +
  geom_line(color = "darkgreen", linewidth = 1.5) +
  geom_point(color = "darkgreen", size = 2) +
  geom_hline(yintercept = K, linetype = "dashed", color = "red") +
  annotate("text", x = 40, y = 108, label = "K = 100", color = "red") +
  labs(x = "Tiempo", y = "N(t)", title = "Crecimiento Logístico Continuo (solución analítica)") +
  theme_bw()

## Comparación: Modelo discreto vs continuo

dlogistic <- function(K = 100, r = 1, N0 = 2, t = 15) {
  N <- c(N0, numeric(t))
  for (i in 1:t) N[i + 1] <- {
    N[i] + r * N[i] * (1 -  N[i]/K)
  }
  return(N)
}

# Modelo continuo (solución analítica)
N_cont <- logistic_solution(t, N0 = 2, r = 1, K = 100)

# Modelo discreto
N_disc <- dlogistic(K = 100, r = 1, N0 = 2, t = 50)

comparison <- data.frame(
  t = 0:50,
  Continuo = N_cont,
  Discreto = N_disc
)

ggplot(comparison, aes(x = t)) +
  geom_line(aes(y = Continuo, color = "Continuo"), linewidth = 1.2) +
  geom_point(aes(y = Discreto, color = "Discreto"), size = 2) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "gray50") +
  labs(x = "Tiempo", y = "N", title = "Comparación: Modelo discreto vs continuo") +
  scale_color_manual(values = c("Continuo" = "blue", "Discreto" = "red")) +
  theme_bw()

## Nota: Para valores pequeños de r (< 0.1), ambos modelos producen resultados muy similares.
## Para r > 1, las diferencias se hacen más notorias.





## Anexo

## Diagramas de Bifurcación

## Un diagrama de bifurcación muestra los valores posibles a largo plazo de N 
## en función del parámetro r. Es una herramienta poderosa para entender 
## la complejidad dinámica de los modelos poblacionales.

num.r <- 201; t <- 400
r.s <- seq(0.5, 2.5, length = num.r)

tmp <- sapply(r.s, function(r) dlogistic(r = r, N0 = 50, K = 100, t = t))

tmp.s <- data.frame(
  N = as.vector(tmp),
  r = rep(r.s, each = t + 1),
  time = rep(0:t, num.r)
)

N.bif <- subset(tmp.s, time > 0.5 * t)

ggplot(N.bif, aes(x = r, y = N)) +
  geom_point(alpha = 0.3, size = 0.5, color = "brown") +
  labs(x = "r (tasa de crecimiento)", y = "N (largo plazo)", title = "Diagrama de Bifurcación - Modelo Logístico Discreto") +
  theme_bw() +
  annotate("text", x = 1.0, y = 5, label = "Un punto\n(estable)", color = "blue", size = 3) +
  annotate("text", x = 1.8, y = 5, label = "Ciclo 2\n(oscilación)", color = "green", size = 3) +
  annotate("text", x = 2.2, y = 5, label = "Ciclo 4", color = "orange", size = 3) +
  annotate("text", x = 2.4, y = 5, label = "Caos", color = "red", size = 3)

## Interpretación:
## - r < 1: Un punto de equilibrio (estable)
## - 1 < r < ~1.7: Oscilaciones damped (convergencia a un punto)
## - ~1.7 < r < ~2.0: Ciclos límite (oscilaciones regulares de período 2, 4, 8...)
## - r > ~2.4: Caos (impredecible a largo plazo)

## Diagramas de tela de araña

## Una forma de visualizar la estabilidad de un modelo poblacional es a través 
## de un diagrama de tela de araña, el cual considera el tamaño actual de la 
## población en el eje X y el tamaño en el tiempo siguiente en el eje Y.

# Ejemplo con r alto (caótico)
out <- dlogistic(N0 = 50, r = 2.5, K = 100, t = 200)
n <- out[101:200]
ntp1 <- n[-1]
nt <- n[-length(n)]

df_cobweb <- data.frame(nt = nt, ntp1 = ntp1)

ggplot(df_cobweb, aes(x = nt, y = ntp1)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_path(color = "gray50", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(y = "N[t+1]", x = "N[t]", title = "Diagrama de tela de araña (r = 2.5, caótico)") +
  theme_bw()

# Comparación con r bajo (estable)
out2 <- dlogistic(N0 = 50, r = 0.5, K = 100, t = 200)
n2 <- out2[101:200]
ntp1_2 <- n2[-1]
nt_2 <- n2[-length(n2)]

df_cobweb2 <- data.frame(nt = nt_2, ntp1 = ntp1_2)

ggplot(df_cobweb2, aes(x = nt, y = ntp1)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_path(color = "gray50", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(y = "N[t+1]", x = "N[t]", title = "Diagrama de tela de araña (r = 0.5, estable)") +
  theme_bw()

## Interpretación de los diagramas de tela de araña:
## - Si la trayectoria converge a la línea diagonal (N[t+1] = N[t]): punto fijo estable
## - Si forma un ciclo (oscila entre valores): dinámica periódica
## - Si llena el espacio de forma aparentemente aleatoria: caos

## Pregunta: ¿Qué diferencia notas entre ambos diagramas?


## Actividad

## Utiliza el set de datos Grizzly contenido en el paquete popbio, calcula los 
## valores de lambda, calcula un lambda promedio y utilizando un modelo 
## logístico con capacidad de carga 100 y como tamaño inicial el tamaño inicial 
## del set de datos Grizzly, genera los gráficos correspondientes del modelo 
## e indica si se genera una dinámica discreta normal, con oscilaciones o caótica. 
## Puedes repetir el proceso utilizando el valor de lambda más grande encontrado.

grizzly <- popbio::grizzly