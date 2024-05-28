#resultados bosques aleatorios
library("RColorBrewer")
library("pheatmap")
library("matrixStats")
library("data.table")
library('randomForest')
library("mlbench")
library("caret")
library("rpart.plot")
library("ggplot2")
library("randomForest")
library("caret")

setwd("C:/Users/Rocio/OneDrive - Universidad de Oviedo/Escritorio/universidad/TFG MATE/Data_TFG_Thyroid_final")
tablesDirectory <- "." 

##################################################################################################################################################
########################################CARGAMOS LOS DATOS############################################################

load("bVals_filtered_12846.rda")
load("bVals_filteredPark_12709.rda")
load("phenoData.rda")
load("phenoDataPark.rda")

#renombramos por facilidad
b<- bVals_filtered
columnas_comunes <- intersect(names(b),names(bTest))
posicionesb<- which(colnames(b) %in% columnas_comunes)
b<- b[,posicionesb]
posicionesbTest<- which(colnames(bTest) %in% columnas_comunes)
bTest<-bTest[,posicionesbTest]

#ahora b y bTest tienen las mismas variables - 12709
btotal<- rbind(b,bTest)

#-----------------------------------------------------------------------------------------------------------------------
#utilizo la muestra filtrada
pvalorestot <- sapply(1:12708,
                      function(i) wilcox.test(btotal[,i]~btotal$Group)$p.value)
names(pvalorestot) <- names(btotal[1:12708])
pvaloresajusttot<-p.adjust(pvalorestot, method="bonferroni", n=12708)
cpgtot<- names(head(sort(pvaloresajusttot),2000))
cpggrouptot<- c(cpgtot,"Group")

b2ktest<- bTest[,cpggrouptot]
b2k<- b[,cpggrouptot]
b2ktotal<- rbind(b2k, b2ktest)
length(intersect(colnames(b2k),colnames(b2ktest)))
#2000 cpg

##################################################################################################################################################
set.seed(10)

#empiezo construyendo un bosque sencillo
start_time <- Sys.time()
rfModel <-randomForest(Group ~ ., data = b2k, nthreads = 12, ntree = 4, cutoff = c(0.7, 0.3), proximity=TRUE, oob.prox=TRUE,keep.tree=TRUE)
end_time <- Sys.time()
tiempo_arbol <- end_time - start_time
print(tiempo_arbol) #tarda minuto y medio

#resultados del bosque:
rfModel$confusion #se equivoca más veces que el ultimo arbol, 8 veces
table(bTest$Group, predict(rfModel, bTest, type="class")) #etiqueta 5 BN mal
print(rfModel) #error OOB 22.86%


#represento los 4 árboles:
reprtree:::plot.getTree(rfModel, k=1)
reprtree:::plot.getTree(rfModel, k=2)
reprtree:::plot.getTree(rfModel, k=3)
reprtree:::plot.getTree(rfModel, k=4)
#no aparecen las variables más importantes de los árboles


#construyo un bosque mas grande
start_time <- Sys.time()
rfModel <-randomForest(Group ~ ., data = b2k, nthreads = 12, ntree = 1000, cutoff = c(0.7, 0.3), proximity=TRUE, oob.prox=TRUE,keep.tree=TRUE,importance=TRUE, margin=TRUE)
end_time <- Sys.time()
tiempo_arbol <- end_time - start_time
print(tiempo_arbol) #3 sec


#resultados:
rfModel$confusion #se equivoca 3 veces, etiqueta 2M mal
t(table(b2ktest$Group, predict(rfModel, b2ktest, type="class"))) #etiqueta 1 M como bn
print(rfModel) #oob 5.13%

varImpPlot(rfModel, type=1)
varImpPlot(rfModel, type=2)


importance_df <- as.data.frame(importance(rfModel))
importance_df <- importance_df[order(-importance_df$MeanDecreaseGini),]
top_variables <- rownames(importance_df)[1:3]
top_variables

#represento algunos árboles
reprtree:::plot.getTree(rfModel, k=4)
reprtree:::plot.getTree(rfModel, k=89)
reprtree:::plot.getTree(rfModel, k=560)


#ahora realizo lo mismo con la muestra test
start_time <- Sys.time()
rfModelT <-randomForest(Group ~ ., data = b2ktest, nthreads = 12, ntree = 1000,  proximity=TRUE, oob.prox=TRUE,keep.tree=TRUE,importance=TRUE, margin=TRUE)
end_time <- Sys.time()
tiempo_arbol <- end_time - start_time
print(tiempo_arbol) #tarda 2.6 sec

rfModelT$confusion #se equivoca 1 vez
t(table(b$Group, predict(rfModelT, b, type="class"))) #etiqueta 3 BN como M
print(rfModelT) #oob 3.57%

rfModelT$importance 
varImpPlot(rfModelT, type=1)#diferentes variables
varImpPlot(rfModelT, type=2)

importance_dfT <- as.data.frame(importance(rfModelT))
importance_dfT <- importance_df[order(-importance_dfT$MeanDecreaseGini),]
top_variablesT <- rownames(importance_dfT)[1:3]
top_variablesT


#construyo un bosque con toda la muestra
start_time <- Sys.time()
rfModelTot<-randomForest(Group ~ ., data = b2ktotal, nthreads = 12, ntree = 1000, cutoff = c(0.7, 0.3), proximity=TRUE, oob.prox=TRUE,keep.tree=TRUE,importance=TRUE, margin=TRUE)
end_time <- Sys.time()
tiempo_arbol <- end_time - start_time
print(tiempo_arbol) #tarda 3 seg

#resultados:
print(rfModelTot) #el error OOB 4.48% 
#los errores de la matriz de confusion -> clasifica 3M como BN
varImpPlot(rfModelTot, type=1)#diferentes variables
varImpPlot(rfModelTot, type=2)


#otros parametros
proximidad <- as.matrix(rfModelTot$proximity)
mds <- cmdscale(1 - proximidad, k = 2)
plot(mds, col = as.numeric(b2ktotal$Group), pch = 16, main = "Gráfico de Proximidad",xlab = "Dimensión 1", ylab = "Dimensión 2")
legend("topright", legend = levels(b2ktotal$Group), col = 1:length(levels(b2ktotal$Group)), pch = 16)

plot(outlier(rfModelTot),type="h", col=c("red","blue")[as.numeric(b2ktotal$Group)], ylab="outlier", main="Medida outliers")
legend("topright", legend = levels(b2ktotal$Group), col = c("red", "blue"), lty = 1) #hay un valor de BN que parece un claro outlier, lo quito?


#represento algunos árboles
reprtree:::plot.getTree(rfModelTot, k=4)
reprtree:::plot.getTree(rfModelTot, k=7)
reprtree:::plot.getTree(rfModelTot, k=500)
reprtree:::plot.getTree(rfModelTot, k=870)

#####################################################################################################################################
#Eliminación recursiva de variables (RFE)

#realizo una vez la eliminación recursiva de variables
rfModel <-randomForest(Group ~ ., data = b2ktotal, ntree = 1000, cutoff = c(0.7, 0.3), importance=TRUE)
varImpPlot(rfModel)
CpGs <- as.data.frame(rfModel$importance)
CpGs$Name <- rownames(CpGs)
CpGs <- CpGs[order(CpGs$MeanDecreaseGini, decreasing = T),]
CpGs <- CpGs[1:10,]
control <- rfeControl(functions=rfFuncs, method="cv", number=67)
bVals_filtered_10 <- btot[,rownames(CpGs)[1:10]] 
bVals_filtered_10$Group <- btot$Group
results <- rfe(b2ktotal[,1:2000], b2ktotal[,2001], rfeControl=control, size=c(1:20))
print(results$variables)
print(results$results)
print(results$optsize) 
print(results$fit)
print(results$optVariables) 


#repito esto 100 veces, para reducir la aleatoriedad
array_acc<- data.frame(Variables=numeric(),
                       Accuracy=numeric(), seed=numeric())
array_imp<- data.frame(BN=numeric(), M=numeric(), MeanDecreaseAccuracy=numeric(), MeanDecreaseGini=numeric(), name=character(), seed=numeric())
matrix_predictors <- data.frame(Name = NA, Seed = NA, Pos = NA) 
semillas <- seq(1,100)

for (i in semillas){
  set.seed(i)
  rfModel <-randomForest(Group ~ ., data = b2ktotal, ntree = 1000, cutoff = c(0.7, 0.3), importance=TRUE)
  #varImpPlot(rfModel)
  CpGs <- as.data.frame(rfModel$importance)
  CpGs$Name <- rownames(CpGs)
  CpGs <- CpGs[order(CpGs$MeanDecreaseGini, decreasing = T),]
  CpGs <- CpGs[1:15,]
  CpGs$seed<-i
  array_imp<-rbind(array_imp, CpGs)
  control <- rfeControl(functions=rfFuncs, method="cv", number=67)
  bVals_filtered_15 <- b2ktotal[,rownames(CpGs)[1:15]] 
  bVals_filtered_15$Group <- b2ktotal$Group
  results <- rfe(bVals_filtered_15[,1:15], bVals_filtered_15[,16], sizes=c(1:15), rfeControl=control)
  #results <- rfe(b2ktotal[,1:2000], b2ktotal[,2001], rfeControl=control, size=c(1:15))
  predict <- c(predictors(results))
  info <- results$results[,1:2]
  info$seed<-i
  array_acc<-rbind(array_acc, info)
  matrix_n <- data.frame(Name = predict, Seed = rep(i, length(predict)), Pos = c(1:length(predict)))
  matrix_predictors <- rbind(matrix_predictors, matrix_n)
  print(plot(results, type=c("g", "o"), main=paste("SEMILLA = ",i)))
}



#precision vs número de variables, siendo la media de las 100 iteraciones
media_accuracy <- aggregate(Accuracy ~ Variables, array_acc, mean)
ggplot(media_accuracy[1:20,], aes(x = Variables, y = Accuracy)) +
  geom_point(color = "blue") +  # Puntos
  geom_line(color = "blue") + # Líneas conectando los puntos
  geom_hline(yintercept = 0.96, linetype = "dashed", color = "red") +
  annotate("text", x = Inf, y = 0.961, label = "Mejor árbol", hjust = 1, vjust = 0) +
  labs(title = "Precisión vs número de variables", x = "Número de variables", y = "Precisión") +
  theme_minimal() +             # Estilo minimalista
  theme(plot.title = element_text(hjust = 0.5))  # Centrar el título del gráfico


#mejores variables segun MDG
length(unique(array_imp$Name)) #85 variables
media_por_name <- aggregate(MeanDecreaseGini ~ Name, data = array_imp, FUN = mean)
array_imp$Name <- factor(array_imp$Name, levels = media_por_name$Name[order(media_por_name$MeanDecreaseGini)])
imp<- subset(array_imp, array_imp$MeanDecreaseGini>0.5)
boxp <- ggplot(imp, aes(x = Name, y = MeanDecreaseGini, fill = Name)) + 
  geom_boxplot() + 
  labs(x = "CpG", y = "MeanDecreaseGini", fill = "Tipo", title = "Mean Decrese Gini para cada variable CpG") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
boxp

#mejores variables segun MDA
impA<- subset(array_imp, array_imp$MeanDecreaseAccuracy>0.006)
boxp2 <- ggplot(impA, aes(x = Name, y = MeanDecreaseAccuracy, fill = Name)) + 
  geom_boxplot() + 
  labs(x = "CpG", y = "MeanDecreaseAccuracy", fill = "Tipo", title = "Mean Decrese Accuracy para cada variable CpG") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
boxp2

#variables mas importantes en la prediccion de M, por MDA
mediaM_por_name <- aggregate(M ~ Name, data = array_imp, FUN = mean)
array_imp$Name <- factor(array_imp$Name, levels = mediaM_por_name$Name[order(mediaM_por_name$M)])
boxp3 <- ggplot(array_imp, aes(x = Name, y = M, fill = Name)) + 
  geom_boxplot() + 
  labs(x = "CpG", y = "M", fill = "Tipo", title = "Importancia de cada CpG al clasificar nódulos malignos")
boxp3

#variables mas importantes en la prediccion de BN, por MDA
mediaBN_por_name <- aggregate(BN ~ Name, data = array_imp, FUN = mean)
array_imp$Name <- factor(array_imp$Name, levels = mediaBN_por_name$Name[order(mediaBN_por_name$BN)])
boxp4 <- ggplot(array_imp, aes(x = Name, y = BN, fill = Name)) + 
  geom_boxplot() + 
  labs(x = "CpG", y = "BN", fill = "Tipo", title = "Importancia de cada Cpg al clasificar nódulos benignos")
boxp4


#extraigo las 13 variables más importantes
posicion_CPG<- aggregate(Pos ~ Name, data = matrix_predictors, FUN = mean)
posicion_CPG <- posicion_CPG[order(posicion_CPG$Pos), ]
cpg <- posicion_CPG[1:13, ]
filtro <- which(posicion_CPG$Name %in% cpg$Name)
mequedocon <- posicion_CPG[filtro,]
datos <- mequedocon[order(mequedocon$Pos), ]
ggplot(datos, aes(x = Pos, y = reorder(Name, - Pos))) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Posición de importancia", y = "Variable CpG", title = "Importancia de las variables CpG")

#construyo bosques con las 13 variables
b13<-b2k[, cpg$Name]
b13$Group<-b2k$Group
btest13<-b2ktest[, cpg$Name]
btest13$Group<-b2ktest$Group

rfModel <-randomForest(Group ~ ., data = btest13, ntree = 1000, cutoff = c(0.7, 0.3), importance=TRUE)
print(rfModel)
rfModel$confusion 
t(table(b13$Group, predict(rfModel, b13, type="class"))) 
print(rfModel) #NO CAMBIAN LOS RESULTADOS


rfModelb <-randomForest(Group ~ ., data = b13, ntree = 1000, cutoff = c(0.7, 0.3), importance=TRUE)
t(table(btest13$Group, predict(rfModelb, btest13, type="class"))) 
print(rfModelb)


rfModeltot <-randomForest(Group ~ ., data = rbind(b13, btest13), ntree = 1000, cutoff = c(0.7, 0.3), importance=TRUE)
print(rfModeltot)
reprtree:::plot.getTree(rfModeltot, k=4)
reprtree:::plot.getTree(rfModeltot, k=89)
reprtree:::plot.getTree(rfModeltot, k=560)

#pruebo a introducir T
b2ktotal$T <- b2ktotal$cg04903825+b2ktotal$cg24227782+b2ktotal$cg26620021
b13$T <- b20$cg04903825+b20$cg24227782+b20$cg26620021
btest13$T <- btest20$cg04903825+btest20$cg24227782+btest20$cg26620021


rfModel2 <-randomForest(Group ~ ., data = btest13, ntree = 1000, cutoff = c(0.7, 0.3), importance=TRUE)
print(rfModel2)
rfModel2$confusion 
t(table(b13$Group, predict(rfModel2, b13, type="class")))
print(rfModel2) #NO CAMBIAN LOS RESULTADOS

rfModelb2 <-randomForest(Group ~ ., data = b13, ntree = 1000, cutoff = c(0.7, 0.3), importance=TRUE)
t(table(btest13$Group, predict(rfModelb2, btest13, type="class"))) #etiqueta 1 m como bn
print(rfModelb2)

#no se utiliza

#compruebo los resultados con 3 variables cpg
btot3<- btotal[, c("cg24227782", "cg04903825" , "cg26620021", "Group")]
btot3$Group <- btotal$Group
rfModel3 <-randomForest(Group ~ ., data = btot3, ntree = 1000, cutoff = c(0.7, 0.3), importance=TRUE)
print(rfModel3)

rfModelb3 <-randomForest(Group ~ ., data = b2k[, c( "cg24227782","cg04903825" , "cg26620021", "Group")], ntree = 1000, cutoff = c(0.7, 0.3), importance=TRUE)
t(table(b2ktest$Group, predict(rfModelb3, b2ktest, type="class")))
print(rfModelb3)

rfModelbp3 <-randomForest(Group ~ ., data = b2ktest[, c("cg24227782", "cg04903825" , "cg26620021","Group")], ntree = 1000, cutoff = c(0.7, 0.3), importance=TRUE)
t(table(b2k$Group, predict(rfModelbp3, b2k, type="class"))) 
print(rfModelbp3)
#igual que con 13 variables