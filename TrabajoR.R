##########################################################################################
#####################                                        			##################
#####################    TRABAJO ORDINARIA BIOINFORMÁTICA R             ##################
#####################                                       			##################
##########################################################################################

### 1.- Carga el archivo de entrada "archivo-entrada.txt" correspondiente. 
### Recuerda incluir los nombres de las columnas como el encabezado.
datos <- read.table("archivo-entrada.txt", header = TRUE, sep = "\t") 
### Examina el archivo con las funciones head(), tail() y dim(). (1 pto)
head (datos) #primeros 6 de la tabla
tail (datos) #ultimos seis de la tabla
dim(datos) #dimensiones que hay 619 filas, 7 columnas

### 2.- En nuestra tabla tenemos los valores de 619 genes para tres réplicas
### del Control y tres réplicas de una Condición. Calcula la media y la desviación
### estándar del control y del tratamiento para cada gen. Puedes emplear cualquier
### función que consideres como apply(). Añade esos valores como columnas nuevas a tu tabla original. 
### Serán 4 columnas nuevas en total, la media de Control, la DE de Control, la media de Condición y la DE de Condición. (1,5 pto)
datos$Media_Control   <- rowMeans(datos[, 2:4])
datos$DE_Control      <- apply(datos[, 2:4], 1, sd)
datos$Media_Condicion <- rowMeans(datos[, 5:7])
datos$DE_Condicion    <- apply(datos[, 5:7], 1, sd)

head(datos)

### 3.- En una sola imagen, muestra un histograma para la media del Control y otro para la media de la Condicíon.
### Pon colores diferentes en cada uno y mantén esos colores para el resto del trabajo. (1 pto).
hist(datos$Media_Control, main = "Histograma media Control",
     col = "tomato", xlab = "Media")
hist(datos$Media_Condicion, main = "Histograma media Condición",
     col = "darkgreen", xlab = "Media")

### 4.- En una sola imagen, muestra un boxplot para la condición y otro para el control. 
### Añade el título, colores y una una leyenda en la esquina superior izquierda. ( 1,5 pto)
boxplot(datos$Media_Control, datos$Media_Condicion,
        col   = c("tomato", "darkgreen"),
        names = c("Control", "Condicion"),
        main  = "Distribución de la media por grupo",
        xlab  = "Grupos", ylab = "Media")

legend("topleft",
       legend = c("Control", "Condicion"),
       fill   = c("tomato", "darkgreen"))


### 5.- Evalúa la normalidad de Control y Condición y acorde con los resultados, haz lo siguiente: (2 pto)
### 		Realiza un test de comparativa de medias para cada gen entre condición y control. 
###			Guarda el p-value de la comparativa y anádelo como una columna extra a tu tabla original.shapiro.test(log2(datos$Media_Control + 1))
shapiro.test(datos$Media_Control)
shapiro.test(datos$Media_Condicion)

datos$pvalue <- sapply(1:nrow(datos), function(i){
wilcox.test(as.numeric(datos [i, 2:4]), as.numeric(datos [i, 5:7]))$p.value
})
head(datos)
table(datos$pvalue)

### 6.- De los resultados del apartado 5, saca (mediante código) un listado de los genes diferencialmente expresados (aquellos que tengan un pvalue < 0.05)
### ¿Cuántos genes tienes? Sino te sale ninguno, justifica por qué. En caso de que sí que te salgan, ¿cuántos salen sobreexpresados y cuántos reprimidos? (1 pto)
deg   <- datos[!is.na(datos$pvalue) & datos$pvalue < 0.05, ]
sobre <- deg[deg$Media_Condicion > deg$Media_Control, ]
repri <- deg[deg$Media_Condicion < deg$Media_Control, ]

cat("Genes diferencialmente expresados:", nrow(deg), "\n")
cat("Sobreexpresados:", nrow(sobre), "\n")
cat("Reprimidos:", nrow(repri), "\n")

deg   <- datos[round(datos$pvalue, 2)<= 0.1, ]
sobre <- deg[deg$Media_Condicion > deg$Media_Control, ]
repri <- deg[deg$Media_Condicion < deg$Media_Control, ]

cat("Genes diferencialmente expresados:", nrow(deg), "\n")
cat("Sobreexpresados:", nrow(sobre), "\n")
cat("Reprimidos:", nrow(repri), "\n")

### 7.- Por último, realiza un heatmap de los genes diferencialmente expresados
### Emplea la Práctica de R como referencia. (2 pts)
mat <- as.matrix(deg[, 2:7])
rownames(mat) <- deg$Gene.Name
colnames(mat) <- c("Control1", "Control2", "Control3",
                   "Condicion1", "Condicion2", "Condicion3")

mat_log <- log2(mat + 1)

# Quitar genes sin variación (darían NaN al escalar)
mat_log <- mat_log[apply(mat_log, 1, sd) > 0, , drop = FALSE]

col_columnas <- rep(c("skyblue", "darkgreen"), each = 3)

heatmap(mat_log,
        scale = "row",
        col = colorRampPalette(c("blue", "white", "red"))(100),
        ColSideColors = col_columnas,
        margins = c(8, 8),
        labRow = if (nrow(mat_log) > 60) NA else rownames(mat_log),
        main = "Heatmap de genes diferencialmente expresados (p <= 0.1)")

legend("topright",
       legend = c("Control", "Condicion"),
       fill = c("skyblue", "darkgreen"),
       bty = "n")

