# Claudia Pérez Marín_Trabajo2.R
# Trabajo final Bioinformática - Curso 25/26
# Análisis de parámetros biomédicos por tratamiento


############################################
# Trabajo Práctico - Análisis Biomed
############################################

# 1. Cargar librerías y datos 
# Instalar si no las tienes:
library(tidyverse)

# Cargar datos (ajusta la ruta si es necesario)
datos <- read.csv("datos_biomed (2).csv", header = TRUE)

# Asegurar que las variables sean del tipo correcto
datos$Glucosa <- as.numeric(datos$Glucosa)
datos$Presion <- as.numeric(datos$Presion)
datos$Colesterol <- as.numeric(datos$Colesterol)
datos$Tratamiento <- as.factor(datos$Tratamiento)

# 2. Exploración inicial --------------------------------------------------
head(datos)
summary(datos)
dim(datos)
str(datos)

num_variables <- ncol(datos)
num_tratamientos <- length(unique(datos$Tratamiento))

cat("Número de variables:", num_variables, "\n")
cat("Número de tratamientos:", num_tratamientos, "\n")

# 3. Boxplot por tratamiento 
ggplot(datos, aes(x = Tratamiento, y = Glucosa, fill = Tratamiento)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Boxplots de Glucosa por Tratamiento")

# 4. Violin plot
ggplot(datos, aes(x = Tratamiento, y = Glucosa, fill = Tratamiento)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "yellow") +
  theme_minimal() +
  labs(title = "Violin plot de Glucosa por Tratamiento")

# 5. Gráfico de dispersión Glucosa vs Presión
plot(datos$Glucosa, datos$Presion,
     col = as.factor(datos$Tratamiento), pch = 19,
     xlab = "Glucosa", ylab = "Presión",
     main = "Dispersión Glucosa vs Presión")

legend("bottomright", legend = levels(datos$Tratamiento),
       col = 1:length(levels(datos$Tratamiento)),
       pch = 19, title = "Tratamiento")

# 6. Facet Grid: Colesterol vs Presión por tratamiento
ggplot(datos, aes(x = Colesterol, y = Presion)) +
  geom_point(aes(color = Tratamiento)) +
  facet_grid(. ~ Tratamiento) +
  theme_minimal() +
  labs(title = "Facet Grid: Colesterol vs Presión por Tratamiento")

# 7. Histogramas para cada variable
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))

hist(datos$Glucosa, main = "Histograma Glucosa", col = "skyblue")
hist(datos$Presion, main = "Histograma Presión", col = "orange")
hist(datos$Colesterol, main = "Histograma Colesterol", col = "green")

par(mfrow = c(1,1))

# 8. Crear un factor a partir del tratamiento
factor_tratamiento <- factor(datos$Tratamiento)
str(factor_tratamiento)

# 9. Media y desviación estándar de glucosa por tratamiento
aggregate(Glucosa ~ Tratamiento, data = datos, FUN = mean)
aggregate(Glucosa ~ Tratamiento, data = datos, FUN = sd)

# 10. Extraer datos de cada tratamiento
lista_tratamientos <- split(datos, datos$Tratamiento)
# Ejemplo: Placebo
placebo <- lista_tratamientos[["Placebo"]]

# 11. Prueba de normalidad y comparativa de medias 
# Normalidad de glucosa por tratamiento
by(datos$Glucosa, datos$Tratamiento, shapiro.test)

# Comparación de ejemplo: Placebo vs FarmacoA (ajusta si es necesario)
if(all(c("Placebo","FarmacoA") %in% levels(datos$Tratamiento))){
  t.test(lista_tratamientos[["Placebo"]]$Glucosa,
         lista_tratamientos[["FarmacoA"]]$Glucosa,
         var.equal = TRUE)
}

# Si los datos no son normales,prueba no paramétrica:
# kruskal.test(Glucosa ~ Tratamiento, data = datos)

# 12. ANOVA sobre glucosa por tratamiento
anova_glucosa <- aov(Glucosa ~ Tratamiento, data = datos)
summary(anova_glucosa)

# Post-hoc Tukey si hay diferencias significativas
TukeyHSD(anova_glucosa)
