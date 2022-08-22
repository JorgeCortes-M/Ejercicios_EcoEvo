library("ggplot2")
library("dplyr")


## Prueba de deriva genética


alelos = c(rep("a",100), rep("A",100))


#Muestreamos una vez

muestra <- sample(alelos, 20, replace = TRUE)

muestra

#Calculo de frecuencias

numA <- sum(muestra=="A")
numa <- sum(muestra=="a")

#numA <- conteo estudiantes
#numa <- conteo estudiantes

num_al <- cbind(numA,numa)
num_al <- as.data.frame(num_al)
num_freqs <- num_al %>% 
  mutate(
    A_nueva = numA/20,
    a_nueva = numa/20,
    a_inicial = 0.5,
    A_inicial = 0.5,
    numA = NULL,
    numa = NULL
  )

#Ordenar los datos

num_freqs <- as.data.frame(num_freqs)
num_freqs<-t(num_freqs)
num_freqs <- as.data.frame(num_freqs)

num_freqs <- num_freqs %>% mutate(
  freq = rownames(num_freqs)
)

#Graficar para comparar frecuencias iniciales y finales

ggplot(data = num_freqs, aes(x=freq, y=V1, fill = freq))+
  geom_bar(position = "identity", stat = "identity", alpha = .5 )+
  ggtitle("Comparación de frecuencias luego de muestreo")+
  xlab("Alelos")+
  ylab("Frecuencia alelos")+
  theme_minimal()


#Muestreo repetidas veces

all_freqA = c()

for(rep in 1:1000){
  nuev_muestra <- sample(alelos, 20, replace = T)
  numA <- sum(nuev_muestra=="A")
  freqA <- numA/20
  all_freqA[rep] <- freqA
}

all_freqA <- as.data.frame(all_freqA)


ggplot(data=all_freqA, aes(x=all_freqA))+
  geom_histogram(binwidth=0.05, colour="black", fill="palegreen2")+
  ggtitle("Frecuencia alelo A en 1000 muestreos")+
  xlab("Frecuencia alelo A")+
  ylab("Conteo")+
  theme_minimal()
