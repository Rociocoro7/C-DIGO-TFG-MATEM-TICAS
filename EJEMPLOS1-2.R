#EJEMPLOS TFG CAPÍTULOS 1-2

library(rpart.plot)
library(rpart)
library(caret)
library(ggplot2)
set.seed(1)

#utilizo la base de datos de iris
summary(iris)

###################################################################################
#EJEMPLOS ÁRBOLES DE CLASIFICACIÓN - datos Iris

#-------------------------GENERALIDADES-------------------
#construyo el primer árbol con el ancho y largo de los pétalos
arbol<-rpart(Species~ Petal.Length+Petal.Width, iris, control=rpart.control(xval=10))
arbol
rpart.plot(arbol, main = "Árbol de clasificación para los datos de Iris", extra=2, faclen=0, type=1, yes.text="si",no.text="no")


#magnitudes del arbol
arbol$frame[,1:5]#los nodos
arbol$where[c(1:5,145:150)]#hojas
table(arbol$where) #tamaño hojas
summary(arbol) #porcentajes

#calcular xerror
arbol$cptable
printcp(arbol)

#matriz de confusion
pred <- predict(arbol, iris, type="class")
confusionMatrix(pred, iris$Species)

#modificamos las variables de entrada
arbol2<- rpart(Species ~ Sepal.Length + Sepal.Width, data=iris)
rpart.plot(arbol2, main = "Árbol de clasificación para los datos de Iris", extra=2, faclen=0, type=1, yes.text="si",no.text="no")
summary(arbol2)
table(iris$Species, predict(arbol,iris,type="class"))


#---------------------------------IMPUREZA--------------------------------------------
#Diferencia entre Gini e información
a.gini<-rpart(Species~ ., iris, parms=list(split="gini"))
prp(a.gini, extra=2, faclen=0, type=0, yes.text="si",no.text="no")
a.info<-rpart(Species~ ., iris, parms=list(split="AUC"))
prp(a.info, extra=2, faclen=0, type=0, yes.text="si",no.text="no")

#representaciones para el ejemplo 2.2.1
rpart.plot(rpart(Species~., iris, control=rpart.control(minsplit = 120)), extra=2, faclen=0, type=1, yes.text="si",no.text="no")
tl<-iris[iris$Petal.Length<4,]
tr<-iris[iris$Petal.Length>=4,]
summary(tl$Species)
summary(tr$Species)

#importancia de las variables
a.gini$variable.importance
a.info$variable.importance
df <- data.frame(importancia = a.gini$variable.importance, "variable" = c("Petal.Width", "Petal.Length", "Sepal.Length","Sepal.Width "))
ggplot2::ggplot(df) + labs(tilte="Importancia de las variables") +
  geom_segment(aes(x = variable, y = importancia, xend = variable, yend = importancia), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = variable, y = importancia, col = variable), 
             size = 4, show.legend = F) +
  coord_flip() +
  theme_bw() + ggtitle("Importancia de las variables")



#---------------------------criterios de parada-------------------------
a1<-rpart(Species~Petal.Length+Petal.Width,iris, control=rpart.control(minsplit=5))
a2<-rpart(Species~Petal.Length+Petal.Width,iris, control=rpart.control(minsplit=55))
prp(a1, extra=2, faclen=0, type=0, yes.text="si",no.text="no")
prp(a2, extra=2, faclen=0, type=0, yes.text="si",no.text="no")

b1<-rpart(Species~Petal.Length+Petal.Width,iris, control=rpart.control(maxdepth=1))
b2<-rpart(Species~Petal.Length+Petal.Width,iris, control=rpart.control(maxdepth=5))
prp(b1, extra=2, faclen=0, type=0, yes.text="si",no.text="no")
prp(b2, extra=2, faclen=0, type=0, yes.text="si",no.text="no")

c1<-rpart(Species~Petal.Length+Petal.Width,iris, control=rpart.control(minbucket=50))
c2<-rpart(Species~Petal.Length+Petal.Width,iris, control=rpart.control(minbucket=5))
prp(c1, extra=2, faclen=0, type=0, yes.text="si",no.text="no")
prp(c2, extra=2, faclen=0, type=0, yes.text="si",no.text="no")


#parámetro de complejidad
printcp(arbol)
plotcp(arbol)

cp_table <- as.data.frame(arbol$cptable)
ggplot(cp_table, aes(x=nsplit, y=xerror)) +
  geom_point(size=3, color="blue") +
  geom_line(aes(group=1), color="blue") +
  geom_errorbar(aes(ymin=xerror - xstd, ymax=xerror + xstd), width=0.2) +
  labs(title="Gráfico de Costos de Complejidad del Árbol de Clasificación",
       x="Número de Divisiones",
       y="Error de Validación Cruzada") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10))

arbolcp<-rpart(Species~.,iris, control=rpart.control(cp=0.001,minsplit=2))
printcp(arbolcp)
rpart.plot(arbolcp, extra=2, faclen=0, type=1, yes.text="si",no.text="no")
plotcp(arbolcp)

cp_table2 <- as.data.frame(arbolcp$cptable)
ggplot(cp_table2, aes(x=nsplit, y=xerror)) +
  geom_point(size=3, color="blue") +
  geom_line(aes(group=1), color="blue") +
  geom_errorbar(aes(ymin=xerror - xstd, ymax=xerror + xstd), width=0.2) +
  labs(title="Gráfico de Costos de Complejidad del Árbol de Clasificación",
       x="Número de Divisiones",
       y="Error de Validación Cruzada") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold"),
        axis.title = element_text(size=12),
        axis.text = element_text(size=10))

######################################################################################
#EJEMPLOS BOSQUE ALEATORIO

library(randomForest)
library(reprtree)

modelo_rf <- randomForest(Species~ ., data = iris, ntree = 3, proximity=TRUE, oob.prox=TRUE, keep.forest=TRUE, importance=TRUE, margin=TRUE)
print(modelo_rf)

reprtree:::plot.getTree(modelo_rf, k=1, text.tree(splits=TRUE, all=TRUE))
reprtree:::plot.getTree(modelo_rf, k=2)
reprtree:::plot.getTree(modelo_rf, k=3)

#hago 10 árboles
modelo_rf10 <- randomForest(Species~ ., data = iris, ntree = 10, proximity=TRUE, oob.prox=TRUE, keep.forest=TRUE, importance=TRUE)
print(modelo_rf10)

#hago 50 árboles
modelo_rf50 <- randomForest(Species~ ., data = iris, ntree = 50, proximity=TRUE, oob.prox=TRUE, keep.forest=TRUE, importance=TRUE)
print(modelo_rf50)

modelo_rf$importance
modelo_rf$err.rate
modelo_rf$confusion
modelo_rf$proximity

#importancia de las variables
importance(modelo_rf,type=1, class=NULL,scale=TRUE)
importance(modelo_rf,type=2,class=NULL,scale=TRUE)

#representacion de la importancia
varImpPlot(modelo_rf, type=1, main="Mean Decrease Accuracy para un bosque de 3 árboles")
varImpPlot(modelo_rf, type=2, main="Mean Decrease Gini para un bosque de 3 árboles")

varUsed(modelo_rf)
getTree(modelo_rf)


rfModel <- randomForest(Species~ ., data = iris, ntree = 1000, proximity=TRUE, oob.prox=TRUE, keep.forest=TRUE, importance=TRUE, margin=TRUE)
print(rfModel)
importance(modelo_rf,type=1, class=NULL,scale=TRUE)
importance(modelo_rf,type=2,class=NULL,scale=TRUE)


plot(margin(rfModel, observed=iris$Species), ylab="Margin")
plot(margin(rfModel), ylab="Margin", pch=19, col="black")

# Añadir puntos coloreados con colores específicos
points(1:length(margins), margins, col=as.numeric(iris$Species), pch=19)
legend("bottomright", legend=levels(iris$Species), col=1:3, pch=19)


predicciones_oob <- predict(modelo_rf, OOB = TRUE)
error_oob <- mean(predicciones_oob != iris$Species)
print(paste("Error OOB:", error_oob))



control <- rfeControl(functions=rfFuncs, method="cv", number=150)
results <- rfe(iris[,1:4], iris[,5], rfeControl=control, size=c(1:4))
results$variables
results$results
results$optsize
results$optVariables
ggplot(results$results, aes(x = Variables, y = Accuracy)) +
  geom_point(color = "blue") +  # Puntos
  geom_line(color = "blue") + # Líneas conectando los puntos
  labs(title = "Precisión vs número de variables", x = "Número de variables", y = "Precisión") +
  theme_minimal() +             # Estilo minimalista
  theme(plot.title = element_text(hjust = 0.5))  # Centrar el título del gráfico
rfModel2 <- randomForest(Species~ Petal.Width + Petal.Length, data = iris, ntree = 1000, proximity=TRUE, oob.prox=TRUE, keep.forest=TRUE, importance=TRUE, margin=TRUE)
print(rfModel2)
outlier_measure <- outlier(rfModel)

# Añadir la medida de outlier al dataframe original
iris$outlier <- outlier_measure

# Plotear los outliers en una gráfica 2D (usando las dos primeras componentes principales)
cmd <- cmdscale(1 - rfModel$proximity, k=2)
cmd_data <- data.frame(cmd, Species=iris$Species, Outlier=iris$outlier)

ggplot(cmd_data, aes(x=V1, y=V2, color=Outlier, shape=Species)) +
  geom_point(size=3) +
  scale_color_gradient(low="blue", high="red") +
  labs(title="Visualización de Outliers usando Random Forest",
       x="Coordenada 1", y="Coordenada 2") +
  theme_minimal() +
  theme(legend.position="right")

plot(iris$outlier)
plot(iris$outlier, main="Medida de Outliers de Iris utilizando un bosque aleatorio",
     xlab="Índice del dato", ylab="Medida de outlier", col=iris$Species,
     pch=19)
legend("topleft", legend=levels(iris$Species), col=1:3, pch=19)

# Instalar y cargar el paquete randomForest si es necesario
if (!requireNamespace("randomForest", quietly = TRUE)) {
  install.packages("randomForest")
}
library(randomForest)

# Usar el dataset iris como ejemplo
data(iris)

# Entrenar un modelo de Random Forest con proximidad
set.seed(123)  # Fijar semilla para reproducibilidad
rf_model <- randomForest(Species ~ ., data=iris, proximity=TRUE)

# Calcular la matriz de proximidad
proximity_matrix <- rfModel$proximity

# Calcular la medida de outlier
n <- nrow(proximity_matrix)
outlier_measure <- numeric(n)

for (i in 1:n) {
  outlier_measure[i] <- sum(proximity_matrix[i, -i]^2)
}

# Añadir la medida de outlier al dataframe original
iris$outlier <- outlier_measure

# Visualizar los outliers con un gráfico simple
plot(iris$outlier, main="Medida de proximidad en Iris",
     xlab="Índice del dato", ylab="Medida de proximidad", col=as.numeric(iris$Species),
     pch=19)
legend("bottomleft", legend=levels(iris$Species), col=1:3, pch=19)
# Determinar las coordenadas para centrar la leyenda
legend_x <- mean(par("usr")[1:2])
legend_y <- mean(par("usr")[3:4])

# Añadir la leyenda en el centro
legend(legend_x, legend_y, legend=levels(iris$Species), col=1:3, pch=19,
       xjust=0.5, yjust=0.5, title="Species", bty="n")

