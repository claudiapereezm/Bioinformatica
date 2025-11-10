#############################################################################
#
# PRACTICA R CLAUDIA PÉREZ MARÍN
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data)
head(data)
tail(data)

# Hacemos un primer histograma para explorar los datos
hist(data, col = "gray", main="GSE5583 - Histogram")

# Transformamos los datos con un logaritmo 
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
# Transformamos nuestro archivo ´data´ a ´data2´para aplicarle la transformación logarítmica con ´log2(data)´. Esto nos permite visualizar mejor los datos, y que se vean con una distribución más homogénea, como ocurre en nuestro caso con el histograma. Aun así, para calcular y hacer los análisis usamos los datos brutos y no la tranformación logarítima, ya que con esta estamos introduciendo errores, por lo que esta última solo la usamos para visualización de gráficos.
data2 = log2(data)
hist(data2, col = "purple", main="GSE5583 (log2) - Histogram")


# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
# ¿Qué es un boxplot?
# Un boxplot (o diagrama de cajas y bigotes) es un tipo de gráfica que permite resumir y visualizar la distribución de un conjunto de datos numéricos. La línea media de los boxes representa la mediana de cada uno, los límites de arriba y abajo de las boxes representan los quartiles Q1 y Q3, los bigotes representan las variaciones hasta los valores máximos y mínimos, y aquellos puntos fuera de los bigotes son los datos que sobresalen del conjunto mayoritario (´outlayers´). Para crear el gráfico, usamos el comando boxplot y lo aplicamos a nuestros datos transformados (data2). Usamos la función col para aplicar color a cada columna en el orden en el que se encuentran en nuestro archivo (por eso salen las 3 columnas de WT en azul y las tres de KO en naranja). Usamos main para representar el título de la gráfica. La función las=2 la usamos para hacer que los nombres en el eje inferior de la gráfica salgan en vertical en lugar de horizontal (se leen de arriba a abajo).
boxplot(data2, col=c("blue", "blue", "blue",
	"orange", "orange", "orange"),
	main="GSE5583 - boxplots", las=2)
	
# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación ç
# de los valores de expresión. ¿Es correcta la separación?
# Si. Estamos clasificando las muestras por patrones de expresión, y usamos este clustering para ver que el experimento y los datos se organizan correctamente. Vemos que los WT y los KO están agrupados juntos entre ellos pero saparados entre cada condición. Usamos hcclust para crearlo y después usamos un comando para la correlación de los datos. Finalmente, hacemos la imagen con el plot, y se imprime el clustering (una representación parecida a un árbol filogenético).
hc = hclust(as.dist(1-cor(data2)))
plot(hc, main="GSE5583 - Hierarchical Clustering")


#######################################
# Análisis de Expresión Diferencial 
#######################################
head(data)
# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
# Se generan dos tablas, una para WT y otra para KO. Para ello, seleccionamos las columnas del archivo orginial (data) que están representadas por las muestras de cada condición. Después usamos `head´para que nos imprima cada tabla.
wt <- data[,1:3]
ko <- data[,4:6]
class(wt)
head(wt)
head(ko)
# Calcula las medias de las muestras para cada condición. Usa apply
wt.mean = apply(wt, 1, mean)
ko.mean = apply(ko, 1, mean)
head(wt.mean)
head(ko.mean)

# ¿Cuál es la media más alta?
# La media más alta es de 37460.5 (calculada con la función ´máx´).
limit = max(wt.mean, ko.mean)
limit

# Ahora hacemos un scatter plot (gráfico de dispersión)
plot(ko.mean ~ wt.mean, xlab = "WT", ylab = "KO",
	main = "GSE5583 - Scatter", xlim = c(0, limit), ylim = c(0, limit))
# Añadir una línea diagonal con abline
abline(0, 1, col = "red")

# ¿Eres capaz de añadirle un grid?
grid()
abline(a, b): línea de pendiente b y ordenada en el origen a
#abline(h=y): línea horizontal
#abline(v=x): línea vertical
abline(1, 2, col = "red")     # línea y = 2x + 1
abline(h = 2, col = "green")  # línea y = 2
abline(v = 3, col = "violet") # línea x = 3

# Calculamos la diferencia entre las medias de las condiciones
diff.mean = wt.mean - ko.mean

# Hacemos un histograma de las diferencias de medias
hist(diff.mean, col = "gray")

# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
# ¿Cuántas valores tiene cada muestra?
# Lo que queremos es hacer los análisis con los datos brutos para evitar posibles sesgos que puedas introducir los datos transformados. Estos solo los vamos a usar al final para visualizar las representaciones gráficas. 
# Cada muestra tiene 12488 valores. Tenemos dos grupos divididos en 3 réplicas cada uno (total de 6 muestras).
pvalue = NULL 
tstat = NULL 
for(i in 1 : nrow(data)) { #Para cada gen
	x = wt[i,] # gene wt número i
	y = ko[i,] # gene ko número i
	
	# Hacemos el test
	t = t.test(x, y)
	
	# Añadimos el p-value a la lista
	pvalue[i] = t$p.value
	# Añadimos las estadísticas a la lista
	tstat[i] = t$statistic
}

head(pvalue)

# Ahora comprobamos que hemos hecho TODOS los cálculos
length(pvalue)

# Hacemos un histograma de los p-values.
# ¿Qué pasa si le ponemos con una transformación de -log10?
# Si hacemos la tranformación la estructura del gráfico y la distribución de los datos cambian, haciendo que la gráfica sea más visual y fácil de comprender. Por ejemplo, en este caso cambia los valores marcados en los ejes de la gráfica, por tanto modificando el grosor y altura de las barras.
hist(pvalue,col="gray")
hist(-log10(pvalue), col = "gray")

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano")

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
diff.mean_cutoff = 2
pvalue_cutoff = 0.01
abline(v = diff.mean_cutoff, col = "blue", lwd = 3)
#abline(v = -diff.mean_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean = abs(diff.mean) >= diff.mean_cutoff
dim(data[filter_by_diff.mean, ])

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue, ])

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
# 426 genes pasan los dos cortes (es decir, están diferencialmente expresados).
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered = data[filter_combined,]
dim(filtered)
head(filtered)

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #2")
points (diff.mean[filter_combined], -log10(pvalue[filter_combined]),col = "red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?
# Los genes sobreexpresados tendrán mayor valor en el KO que en el WT, por lo que tendrán signos negativos, y por eso salen los sobreexpresados en el lado izquierdo de la gráfica. Esto es porque para calcular la diferencia de medias usamos la resta WT - KO.
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #3")
points (diff.mean[filter_combined & diff.mean < 0],
	-log10(pvalue[filter_combined & diff.mean < 0]), col = "red")
points (diff.mean[filter_combined & diff.mean > 0],
	-log10(pvalue[filter_combined & diff.mean > 0]), col = "blue")


# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
# La función heatmap (mapa de calor) hace una representación visual de una matriz de datos, como en nuestro caso en expresión génica. En el código también usamos un clustering jerárquico para ordenar filas y columnas. Con ´filtered´ usamos la matriz de datos que se va a representar; con ´row´ y ´colv´definimos los dendogramas que ordenan las filas y columnas mediante clustering jerárquico, calculado a partir de las correlaciones matriz (cor(), as.dist(), hclust(), as.dendogram(); el argumento cexCol=0.7 ajusta el tamaño de las etiquetas de las columnas y labRow=False indica que no se muestren las etiquetas de las filas.
#  Se puede modificar la paleta de colores utilizando el argumento co1 dentro de la función heatmap(). 

rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,labRow=FALSE)

heatmap(filtered)


# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)
library(RColorBrewer)

# Hacemos nuestro heatmap
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row")

# Lo guardamos en un archivo PDF
pdf ("GSE5583_DE_Heatmap.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row",labRow=FALSE)
dev.off()
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,col = redgreen(75), scale = "row",labRow=FALSE)

# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered, "GSE5583_DE.txt", sep = "\t",quote = FALSE)
