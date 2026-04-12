# Ejercicios de Ecología y Evolución

Colección de talleres prácticos en R sobre temas de ecología de poblaciones y evolución biológica.

## Estructura del Proyecto

```
Ejercicios_EcoEvo/
├── Talleres/
│   ├── Taller_1/          # Crecimiento poblacional
│   ├── Taller_2/          # Poblaciones estructuradas
│   ├── Taller_3/          # Estocasticidad ambiental
│   └── Taller_4/          # Interacciones entre especies
├── Ejercicio_deriva.R     # Ejercicio de deriva genética
└── README.md
```

## Contenido

### Taller 1: Crecimiento Poblacional
- **Crecimiento exponencial**: Modelos de crecimiento sin límites
- **Crecimiento logístico**: Modelos con capacidad de carga (K)
- Archivos: `Crecimiento_exponencial.R`, `Crecimiento_Logistico.R`
- Notebooks:formatos `.Rmd` y `.html`

### Taller 2: Poblaciones Estructuradas
- **Matrices de población**: Modelos demográficos con clases de edad/estado
- **Tablas de vida**: Análisis de supervivencia y fecundidad
- **Análisis demográfico**: Datos de *Vulpes* (zorro ártico) y *Gorilla*
- Archivos de datos: `matrix_vulpes_*.csv`, `Gorilla_lt.xlsx`
- Referencia: Bronikowski et al. (2016) *Sci Data*

### Taller 3: Estocasticidad Ambiental
- **Modelos estocásticos**: Simulaciones con variabilidad ambiental
- Archivo: `Estoc_ambiental.R`

### Taller 4: Interacciones
- **Simulaciones interespecíficas**: Competencia y depredación
- Archivo: `Interact_script.R`

### Ejercicios Adicionales
- **Deriva genética**: Efectos del muestreo finito en frecuencias alélicas
- Archivo: `Ejercicio_deriva.R`

## Requisitos

- **R** (versión 3.6+)
- **Paquetes**: `deSolve`, `popbio`, `xlsx`, `tidyverse`, `ggplot2`

## Uso

1. Abrir los scripts `.R` en RStudio o ejecutarlos en consola
2. Los notebooks `.Rmd` pueden abrirse en RStudio para visualización interactiva
3. Los archivos `.html` permiten ver los notebooks renderizados sin necesidad de R
4. Los datos están en formatos `.csv` y `.xlsx` para análisis adicionales