############################################################################
##              POBLACIONES ESTRUCTURADAS: MATRICES DEMOGRAFICAS          ##
############################################################################
##
## Es posible describir la dinamica poblacional de variados organismos al
## clasificar a sus individuos en "estados" (clases de edad, etapas de
## desarrollo, tallas, etc.) y modelar las transiciones entre ellos usando
## matrices de proyeccion.
##
## Este script cubre:
##   1. Tablas de vida y parametros demograficos (R0, T, r, lambda)
##   2. Construccion de la Matriz de Leslie (por edad)
##   3. Proyecciones poblacionales
##   4. Eigenanalisis completo: lambda, w, v, rho
##   5. Sensibilidad y elasticidad
##   6. Matrices de Lefkovitch (por etapa)
##   7. Ejercicio: Vulpes vulpes + LTRE
##
## INSTRUCCION: todos los archivos de datos deben estar en la misma
## carpeta que este script. No es necesario definir setwd().
############################################################################


## ---- 1. PAQUETES --------------------------------------------------------

packages <- c("ggplot2", "popdemo", "readr", "dplyr", "tidyr",
              "igraph", "Rage", "readxl", "knitr")

installed_packages <- packages %in% rownames(installed.packages())
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

library(popdemo)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)
library(Rage)
library(readxl)
library(knitr)


## ---- 2. TABLA DE VIDA ---------------------------------------------------

## Columnas: x (edad), N (abundancia), lx (supervivencia), mx (fecundidad)
life_table <- read_csv("Life_tab_2.csv", show_col_types = FALSE)
life_table


## ---- 3. PARAMETROS DEMOGRAFICOS -----------------------------------------

## R0: tasa neta de reproduccion      = sum(lx * mx)
## T:  tiempo generacional            = sum(x * lx * mx) / R0
## r:  tasa intrinseca de incremento  = ln(R0) / T
## lambda = e^r  (tasa finita de crecimiento)

R0_val <- sum(life_table$lx * life_table$mx)
T_val  <- sum(life_table$x * life_table$lx * life_table$mx) / R0_val
r_val  <- log(R0_val) / T_val

cat("R0     =", round(R0_val, 4), "\n")
cat("T      =", round(T_val,  4), "anos\n")
cat("r      =", round(r_val,  4), "por ano\n")
cat("lambda = e^r =", round(exp(r_val), 4), "\n\n")

## Verificacion de la relacion R0 - lambda - T:
cat("R0 desde lambda y T: lambda^T =", round(exp(r_val)^T_val, 4),
    " | R0 observado =", round(R0_val, 4), "\n")

## Agregar columnas auxiliares a la tabla
life_table <- life_table %>%
  mutate(
    lx_mx   = lx * mx,
    x_lx_mx = x  * lx * mx
  )

life_table


## ---- 4. CURVA DE SUPERVIVENCIA ------------------------------------------

ggplot(life_table, aes(x = x, y = lx)) +
  geom_line(color = "firebrick", linewidth = 1) +
  geom_point(color = "firebrick", size = 2.5) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title = "Curva de supervivencia",
    x     = "Clase de edad (x)",
    y     = "lx"
  ) +
  theme_classic(base_size = 13)


## ---- 5. MATRIZ DE LESLIE ------------------------------------------------

## Modelo de nacimiento pulsado con censo post-nacimiento:
##   Gi = lx[i+1] / lx[i]       probabilidad de supervivencia
##   Fi = Gi * mx[i+1]           fertilidad de la clase i

n  <- nrow(life_table)
lx <- life_table$lx
mx <- life_table$mx

Gi <- lx[2:n] / lx[1:(n - 1)]
Fi <- Gi * mx[2:n]

data.frame(Clase = seq_len(n - 1), Gi = round(Gi, 4), Fi = round(Fi, 4))


## Funcion para construir la Matriz de Leslie a partir de Gi y Fi

leslie_from_table <- function(Gi, Fi) {
  n <- length(Gi)
  A <- matrix(0, nrow = n, ncol = n)
  A[1, ] <- Fi                                      # primera fila: fertilidades
  for (i in seq_len(n - 1)) A[i + 1, i] <- Gi[i]  # subdiagonal: supervivencias
  A
}

A <- leslie_from_table(Gi, Fi)
round(A, 4)


## Grafo del ciclo de vida
plot_life_cycle(A, shape = "circle", edgecol = "firebrick")


## ---- 6. PROYECCIONES POBLACIONALES --------------------------------------

## Funcion de proyeccion reutilizable para cualquier matriz y horizonte temporal

proyectar <- function(A, N0, anos) {
  n <- nrow(A)
  N <- matrix(0, nrow = n, ncol = anos + 1)
  N[, 1] <- N0
  for (t in seq_len(anos)) N[, t + 1] <- A %*% N[, t]
  N
}

## Vector inicial de abundancias (excluye la fila de edad 0)
N0    <- as.numeric(life_table$N[-1])
years <- 10

N_proj <- proyectar(A, N0, years)

## Proyeccion por clase
N_df <- as.data.frame(t(N_proj))
colnames(N_df) <- paste0("Clase_", seq_len(nrow(A)))
N_df$Ano <- 0:years

pivot_longer(N_df, cols = -Ano, names_to = "Clase", values_to = "N") %>%
  ggplot(aes(x = Ano, y = N, color = Clase)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  labs(
    title = "Proyeccion poblacional por clase de edad",
    x     = "Ano",
    y     = "Numero de individuos"
  ) +
  theme_classic(base_size = 13)


## Distribucion estable de edades (proporciones)
N_totals <- colSums(N_proj)
N_prop   <- sweep(N_proj, 2, N_totals, "/")

N_prop_df <- as.data.frame(t(N_prop))
colnames(N_prop_df) <- paste0("Clase_", seq_len(nrow(A)))
N_prop_df$Ano <- 0:years

pivot_longer(N_prop_df, cols = -Ano, names_to = "Clase", values_to = "Proporcion") %>%
  ggplot(aes(x = Ano, y = Proporcion, color = Clase)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  labs(
    title = "Convergencia a la distribucion estable de edades",
    x     = "Ano",
    y     = "Proporcion (N / N_total)"
  ) +
  theme_classic(base_size = 13)


## ---- 7. EIGENANALISIS ---------------------------------------------------

eigs_A  <- eigen(A)
dom_pos <- which.max(Re(eigs_A$values))

## lambda1: eigenvalor dominante = tasa finita de crecimiento
lambda1 <- Re(eigs_A$values[dom_pos])

## w: eigenvector derecho = distribucion estable de edades
w      <- Re(eigs_A$vectors[, dom_pos])
w_norm <- w / sum(w)

## v: eigenvector izquierdo = valor reproductivo relativo
eigs_left <- eigen(t(A))
dom_pos_l <- which.max(Re(eigs_left$values))
v         <- Re(eigs_left$vectors[, dom_pos_l])
v_norm    <- v / v[1]   # normalizar: clase 1 = 1

## rho: ratio de amortiguamiento = velocidad de convergencia al estado estable
lambda2 <- sort(Mod(eigs_A$values), decreasing = TRUE)[2]
rho     <- lambda1 / lambda2

cat("lambda1 =", round(lambda1, 6), "\n")
cat("r       =", round(log(lambda1), 6), "por ano\n")
cat("rho     =", round(rho, 4), "\n\n")

data.frame(
  Clase             = paste0("Clase_", seq_len(nrow(A))),
  w_estable         = round(w_norm, 4),
  v_reproductivo    = round(v_norm, 4)
)


## Distribucion estable de edades
ggplot(data.frame(Clase = paste0("Clase_", seq_along(w_norm)), Proporcion = w_norm),
       aes(x = Clase, y = Proporcion)) +
  geom_col(fill = "steelblue", alpha = 0.85) +
  labs(title = "Distribucion estable de edades (w)",
       x = "Clase de edad", y = "Proporcion") +
  theme_classic(base_size = 13)

## Valor reproductivo
ggplot(data.frame(Clase = paste0("Clase_", seq_along(v_norm)), ValorReprod = v_norm),
       aes(x = Clase, y = ValorReprod)) +
  geom_col(fill = "coral3", alpha = 0.85) +
  labs(title = "Valor reproductivo relativo (v)",
       subtitle = "Relativo a la Clase 1 = 1",
       x = "Clase de edad", y = "Valor reproductivo") +
  theme_classic(base_size = 13)


## Convergencia de lambda_t hacia lambda1 en proyeccion larga

years_conv <- 100
N_conv     <- proyectar(A, N0, years_conv)
N_tot_conv <- colSums(N_conv)
Rs_conv    <- N_tot_conv[-1] / N_tot_conv[-(years_conv + 1)]

ggplot(data.frame(Ano = seq_len(years_conv), Lambda_t = Rs_conv),
       aes(x = Ano, y = Lambda_t)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(color = "steelblue", size = 0.8) +
  geom_hline(yintercept = lambda1, color = "firebrick",
             linetype = "dashed", linewidth = 1) +
  annotate("text", x = 80, y = lambda1 + 0.004,
           label = paste0("lambda1 = ", round(lambda1, 4)),
           color = "firebrick", size = 4) +
  labs(
    title = "Convergencia de la tasa de crecimiento anual hacia lambda1",
    x     = "Ano de proyeccion",
    y     = "lambda_t = N(t+1) / N(t)"
  ) +
  theme_classic(base_size = 13)


## ---- 8. SENSIBILIDAD Y ELASTICIDAD --------------------------------------

## Sensibilidad: cambio absoluto en lambda1 ante cambio absoluto en a_ij
S <- sens(A, eval = "max")

as.data.frame(S) %>%
  setNames(paste0("Clase_", seq_len(ncol(S)))) %>%
  mutate(Origen = paste0("Clase_", seq_len(nrow(S)))) %>%
  pivot_longer(cols = -Origen, names_to = "Destino", values_to = "Sensibilidad") %>%
  ggplot(aes(x = Destino, y = Origen, fill = Sensibilidad)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Sensibilidad de lambda1",
       x = "Destino (columna)", y = "Origen (fila)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## Elasticidad: cambio proporcional en lambda1 ante cambio proporcional en a_ij
## Propiedad clave: los valores de elasticidad siempre suman 1
E <- elas(A, eval = "max")

as.data.frame(E) %>%
  setNames(paste0("Clase_", seq_len(ncol(E)))) %>%
  mutate(Origen = paste0("Clase_", seq_len(nrow(E)))) %>%
  pivot_longer(cols = -Origen, names_to = "Destino", values_to = "Elasticidad") %>%
  ggplot(aes(x = Destino, y = Origen, fill = Elasticidad)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "coral3") +
  labs(title = "Elasticidad de lambda1",
       x = "Destino (columna)", y = "Origen (fila)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cat("Suma de elasticidades:", round(sum(E), 6), "(debe ser = 1)\n")


## ---- 9. MATRIZ DE LEFKOVITCH (por etapa) --------------------------------

## A diferencia de la Matriz de Leslie, la de Lefkovitch permite:
##   - Permanencia en la misma etapa (elementos en la diagonal principal)
##   - Clasificacion por etapa de desarrollo en vez de por edad exacta
## Es la forma mas comun en estudios empiricos de plantas, reptiles y peces.

## Ejemplo: planta perenne con 3 etapas (Plantula, Juvenil, Adulto)

etapas <- c("Plantula", "Juvenil", "Adulto")

B <- matrix(
  c(0.000, 0.000, 4.500,
    0.250, 0.350, 0.000,
    0.000, 0.550, 0.750),
  nrow = 3, ncol = 3, byrow = TRUE,
  dimnames = list(etapas, etapas)
)

## Interpretacion:
##   B[1,3] = 4.5:  cada adulto produce 4.5 plantulas/anio
##   B[2,1] = 0.25: 25% de plantulas maduran a juvenil
##   B[2,2] = 0.35: 35% de juveniles permanecen como juveniles
##   B[3,2] = 0.55: 55% de juveniles maduran a adulto
##   B[3,3] = 0.75: 75% de adultos sobreviven y permanecen como adultos

round(B, 3)

## Grafo del ciclo de vida (nótese el lazo de permanencia)
plot_life_cycle(B, shape = "circle", edgecol = "darkgreen")

## Eigenanalisis de la Matriz de Lefkovitch
eigs_B    <- eigen(B)
dom_B     <- which.max(Re(eigs_B$values))
lambda_B  <- Re(eigs_B$values[dom_B])

w_B       <- Re(eigs_B$vectors[, dom_B])
w_B_norm  <- w_B / sum(w_B)

eigs_B_l  <- eigen(t(B))
v_B       <- Re(eigs_B_l$vectors[, which.max(Re(eigs_B_l$values))])
v_B_norm  <- v_B / v_B[1]

lambda2_B <- sort(Mod(eigs_B$values), decreasing = TRUE)[2]
rho_B     <- lambda_B / lambda2_B

cat("lambda1 (planta perenne):", round(lambda_B, 4), "\n")
cat("r equivalente:           ", round(log(lambda_B), 4), "por ano\n")
cat("rho (amortiguamiento):   ", round(rho_B, 4), "\n\n")

data.frame(
  Etapa          = etapas,
  w_estable      = round(w_B_norm, 4),
  v_reproductivo = round(v_B_norm, 4)
)

## Elasticidad de la Matriz de Lefkovitch
E_B <- elas(B, eval = "max")

as.data.frame(E_B) %>%
  setNames(etapas) %>%
  mutate(Origen = etapas) %>%
  pivot_longer(cols = -Origen, names_to = "Destino", values_to = "Elasticidad") %>%
  ggplot(aes(x = Destino, y = Origen, fill = Elasticidad)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Elasticidad, 3)), size = 4.5) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  labs(title = "Elasticidad (Lefkovitch - planta perenne)",
       x = "Etapa destino", y = "Etapa origen") +
  theme_minimal(base_size = 12)

cat("Suma de elasticidades:", round(sum(E_B), 6))


## ---- 10. EJERCICIO: VULPES VULPES ---------------------------------------

## Dos matrices para Vulpes vulpes en Australia (COMADRE).
## Fuente: Devenish-Nelson et al. (2013), Oikos. IDs: 249920 y 249921.

clases_v <- c("0+", "1+", "2+", ">=3")

vulpes_haunt <- matrix(
  c(0.370, 0.610, 1.210, 1.580,
    0.300, 0.000, 0.000, 0.000,
    0.000, 0.350, 0.000, 0.000,
    0.000, 0.000, 0.570, 0.700),
  nrow = 4, ncol = 4, byrow = TRUE,
  dimnames = list(clases_v, clases_v)
)

vulpes_unman <- matrix(
  c(0.686, 1.271, 1.426, 0.332,
    0.390, 0.000, 0.000, 0.000,
    0.000, 0.650, 0.000, 0.000,
    0.000, 0.000, 0.920, 0.180),
  nrow = 4, ncol = 4, byrow = TRUE,
  dimnames = list(clases_v, clases_v)
)

## Funcion auxiliar: eigenanalisis completo en un solo objeto
calc_eigen_completo <- function(mat) {
  eig   <- eigen(mat)
  dpos  <- which.max(Re(eig$values))
  lam   <- Re(eig$values[dpos])
  w_raw <- Re(eig$vectors[, dpos])
  eig_l <- eigen(t(mat))
  v_raw <- Re(eig_l$vectors[, which.max(Re(eig_l$values))])
  list(
    lambda = lam,
    w      = w_raw / sum(w_raw),
    v      = v_raw / v_raw[1],
    rho    = lam / sort(Mod(eig$values), decreasing = TRUE)[2]
  )
}

res_haunt <- calc_eigen_completo(vulpes_haunt)
res_unman <- calc_eigen_completo(vulpes_unman)

## Comparacion de parametros demograficos
data.frame(
  Parametro = c("lambda1", "r", "rho"),
  Caza      = c(round(res_haunt$lambda, 4),
                round(log(res_haunt$lambda), 4),
                round(res_haunt$rho, 4)),
  Control   = c(round(res_unman$lambda, 4),
                round(log(res_unman$lambda), 4),
                round(res_unman$rho, 4))
)

## Distribucion estable y valor reproductivo
data.frame(
  Clase   = clases_v,
  w_haunt = round(res_haunt$w, 4),
  w_unman = round(res_unman$w, 4),
  v_haunt = round(res_haunt$v, 4),
  v_unman = round(res_unman$v, 4)
)

## Grafos del ciclo de vida
par(mfrow = c(1, 2))
plot_life_cycle(vulpes_haunt, shape = "circle", edgecol = "firebrick")
title(main = "Caza intensa", cex.main = 1.2)
plot_life_cycle(vulpes_unman, shape = "circle", edgecol = "steelblue")
title(main = "Control", cex.main = 1.2)
par(mfrow = c(1, 1))

## Proyeccion a 20 anos
N0_v    <- c(100, 50, 20, 10)
years_v <- 20

Nh <- proyectar(vulpes_haunt, N0_v, years_v)
Nu <- proyectar(vulpes_unman, N0_v, years_v)

data.frame(
  Ano          = 0:years_v,
  Caza_intensa = colSums(Nh),
  Control      = colSums(Nu)
) %>%
  pivot_longer(cols = -Ano, names_to = "Tratamiento", values_to = "N_total") %>%
  ggplot(aes(x = Ano, y = N_total, color = Tratamiento)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Caza_intensa" = "firebrick", "Control" = "steelblue")) +
  labs(title = "Vulpes vulpes - Proyeccion a 20 anos",
       x = "Ano", y = "Abundancia total (N)") +
  theme_classic(base_size = 13)

## Elasticidades comparativas
plot_elas <- function(E, titulo, clases) {
  as.data.frame(E) %>%
    setNames(clases) %>%
    mutate(Origen = clases) %>%
    pivot_longer(cols = -Origen, names_to = "Destino", values_to = "Elasticidad") %>%
    ggplot(aes(x = Destino, y = Origen, fill = Elasticidad)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Elasticidad, 3)), size = 4) +
    scale_fill_gradient(low = "white", high = "coral3") +
    labs(title = titulo, x = "Clase destino", y = "Clase origen") +
    theme_minimal(base_size = 11)
}

E_haunt <- elas(vulpes_haunt, eval = "max")
E_unman <- elas(vulpes_unman, eval = "max")

print(plot_elas(E_haunt, "Elasticidad - Caza intensa", clases_v))
print(plot_elas(E_unman, "Elasticidad - Control",      clases_v))


## ---- 11. LTRE (Life Table Response Experiment) --------------------------

## Descompone la diferencia en lambda entre dos poblaciones en contribuciones
## por elemento de la matriz.
## Referencia: matriz promedio de las dos poblaciones.
## Contribucion de a_ij: c_ij = (a_ij_haunt - a_ij_unman) * s_ij_referencia

A_ref  <- (vulpes_haunt + vulpes_unman) / 2
S_ref  <- sens(A_ref, eval = "max")
C_ltre <- (vulpes_haunt - vulpes_unman) * S_ref

cat("Diferencia observada en lambda (haunt - unman):",
    round(res_haunt$lambda - res_unman$lambda, 6), "\n")
cat("Suma de contribuciones LTRE:                  ",
    round(sum(C_ltre), 6), "\n")

as.data.frame(C_ltre) %>%
  setNames(clases_v) %>%
  mutate(Origen = clases_v) %>%
  pivot_longer(cols = -Origen, names_to = "Destino", values_to = "Contribucion") %>%
  ggplot(aes(x = Destino, y = Origen, fill = Contribucion)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Contribucion, 3)), size = 4) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick", midpoint = 0) +
  labs(
    title    = "Contribuciones LTRE a Delta-lambda (haunt - unman)",
    subtitle = "Rojo: la caza aumenta lambda | Azul: la caza disminuye lambda",
    x        = "Clase destino",
    y        = "Clase origen"
  ) +
  theme_minimal(base_size = 12)


## ---- 12. BONUS: GORILLA GORILLA -----------------------------------------

## Tabla de vida para Gorilla gorilla.
## Tarea: calcular R0, T, r; construir Matriz de Leslie; eigenanalisis completo;
## sensibilidad y elasticidad. Comparar con los resultados de Vulpes vulpes.

Gorilla <- read_xlsx("Gorilla_lt.xlsx", sheet = "Hoja1")
Gorilla

## Su analisis aqui:


## ---- SESSION INFO -------------------------------------------------------
sessionInfo()
