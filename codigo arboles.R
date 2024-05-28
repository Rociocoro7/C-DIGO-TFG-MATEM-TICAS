library("RColorBrewer")
library("pheatmap")
library("matrixStats")
library("data.table")
library('randomForest')
library("mlbench")
library("caret")
library("rpart.plot")
library("ggplot2")

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
bTest<- bVals_filteredPark

#Nos quedamos con las CpG comunes a ambas muestras
columnas_comunes <- intersect(names(b),names(bTest))
posicionesb<- which(colnames(b) %in% columnas_comunes)
b<- b[,posicionesb]
posicionesbTest<- which(colnames(bTest) %in% columnas_comunes)
bTest<-bTest[,posicionesbTest]

#ahora b y bTest tienen las mismas variables
btotal<- rbind(b,bTest)

##################################################################################################################################################
set.seed(10)

#muestras de datos:
summary(b$Group) #17 M y 22 BN -> 0.44 0.56%
summary(bTest$Group) #7BN 21 M -> 0.25 y 0.75%
summary(btotal$Group) #29 BN 38M 

summary(b$cg17154646)
summary(b$cg16177576)

by(b$cg17154646 , b$Group , mean)
by(b$cg17154646 , b$Group , sd)
by(b$cg17154646 , b$Group , min)
by(b$cg17154646 , b$Group , max)


#representamos el valor de las CpG para una única biopsia, en particular, la biopsia 24
ggplot(data = data.frame(x = bValues[,"BIOP_24"]), aes(x = x)) +
  geom_histogram(binwidth = 0.01, fill = "lightblue", color = "lightblue", aes(y=..density..), alpha=0.8) +
  theme_minimal() +
  labs(title = "Distribución del valor de las CpG para la biopsia 24",
       x = "Valores de las CpG",
       y = "Frecuencia") +
  theme(plot.title = element_text(hjust = 0.5), # Centrar el título
        axis.text = element_text(size = 10),    # Tamaño del texto en los ejes
        axis.title = element_text(size = 12),   # Tamaño del texto de los títulos de los ejes
        legend.position = "none")

#representamos dos CpG, una potencialmente discriminante y otra no

boxp = ggplot(btotal, aes(x = factor(Group), y =cg17154646 , fill = Group)) + 
  geom_boxplot()+ labs(x = "Tipo Tumor", y = "cg17154646",fill='Tipo', title = "CpG 17154646")
boxp

boxp2 = ggplot(btotal, aes(x = factor(Group), y =cg16177576 , fill = Group)) + 
  geom_boxplot()+ labs(x = "Tipo Tumor", y = "cg16177576",fill='Tipo', title = "CpG 16177576")
boxp2



#represento la distribución de una CpG para todas las biopsias, en particular la CpG24227782

ggplot(data = data.frame(x = btotal[,"cg24227782"]), aes(x = x)) +
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "lightblue", aes(y=..density..), alpha=0.8) +
  theme_minimal() +
  labs(title = "Distribución del valor de la CpG24227782",
       x = "CpG24227782",
       y = "Frecuencia") +
  theme(plot.title = element_text(hjust = 0.5), # Centrar el título
        axis.text = element_text(size = 10),    # Tamaño del texto en los ejes
        axis.title = element_text(size = 12),   # Tamaño del texto de los títulos de los ejes
        legend.position = "none")

###########################################################################
#definimos funciones que posteriormente se utilizarán

#Función para calcular que variable de las más importantes clasifica mejor la muestra test
evaluar_variables <- function(variables, datos_entrenamiento, datos_prueba) {
  error <- c()
  fallos <- c()
  n<- nrow(datos_entrenamiento)
  for (var in variables) {
    formula <- as.formula(paste("Group ~", var))
    arbol_var <- rpart(formula, datos_entrenamiento, control = rpart.control(xval = n, cp = 0))
    pred <- predict(arbol_var, datos_prueba, type = 'class')
    t <- table(datos_prueba$Group, pred)
    precision <- sum(diag(t)) / sum(t)
    fallosi <- sum(t) - sum(diag(t))
    
    error <- c(error, precision)
    fallos <- c(fallos, fallosi)
  }
  
  resumen <- data.frame(variables = variables, precision = error, fallos = fallos)
  return(resumen)
}

###############################################################################
##################################################################################
#################################################################################
############## ÁRBOLES DE CLASIFICACIÓN ########################################

#empezamos con un árbol de clasificación con 39 xval, entrenado con la muestra bValues
#en cada validación queda 1 muestra fuera 
#cp=0, sin limite de profundidad

start_time <- Sys.time()
arbol<-rpart(Group~.,b,control=rpart.control(xval=39, cp=0))
end_time <- Sys.time()
tiempo_arbol <- end_time - start_time
print(tiempo_arbol)
#tarda 34.5314s

#resultados: 
print(arbol) #cg22956858
printcp(arbol)
plotcp(arbol)
#xerror 0.4359*0.11765=0.05
rpart.plot(arbol, main = "Árbol de clasificación con la muestra bValues", extra=2, faclen=0, type=1, yes.text="si",no.text="no")
#sólo utiliza cg22956858
#aciertos y fallos del árbol
pred<-predict(arbol,b, type='class')
t<- table(b$Group, pred)
prop.table(t,1) #probabilidad de fallar cada tipo
t(t) #clasificamos 1 BN como M
#las filas son las clases predichas 

#evaluamos en la muestra test
predTest<-predict(arbol,bTest, type='class')
tTest<- table(bTest$Group, predTest)
t(tTest) #9 fallos, clasifica 9 malignos como benignos
#el error 9/28=0.32

arbol$variable.importance #hay tres cpg igual de importantes cg04903825 cg05372993 cg22956858
cor(b[,c("cg04903825","cg05372993", "cg22956858")])
#son variables muy correlacionadas, les siguen otras tres cpg igua de importantes cg15625618 cg20089319 cg21137882 
df <- data.frame(importancia = arbol$variable.importance, "variable" = c("cg04903825", "cg05372993", "cg22956858","cg15625618","cg20089319","cg21137882"))
ggplot2::ggplot(df) + labs(tilte="Importancia de las variables") +
  geom_segment(aes(x = variable, y = importancia, xend = variable, yend = importancia), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = variable, y = importancia, col = variable), 
             size = 4, show.legend = F) +
  coord_flip() +
  theme_bw() + ggtitle("Importancia de las variables")

variables <- names(arbol$variable.importance)
evaluar_variables(variables, b, bTest)
#La mejor es cg04903825

arbol2<-rpart(Group~cg04903825,b,control=rpart.control(xval=39, cp=0))
predTest<-predict(arbol2,bTest, type='class')
tTest<- table(bTest$Group, predTest)
t(tTest) 

#repetimos el proceso para la muestra test
start_time<- Sys.time()
arbolTest<-rpart(Group~.,bTest,control=rpart.control(xval=28, cp=0))
end_time <- Sys.time()
tiempo_arbol <- end_time - start_time
print(tiempo_arbol) #33.5228 s

#resultados:
print(arbolTest) #cg00175838
rpart.plot(arbolTest, main = "", extra=2, faclen=0, type=1, yes.text="si",no.text="no")
printcp(arbolTest)
#xerror 0.25*0.57=0.14

arbolTest$variable.importance #hay 5 cpg igual de importantes cg00175838 cg01890836 cg04596071 cg13235366 cg13318914 cg23945952
cor(bTest[,c("cg00175838", "cg01890836", "cg04596071", "cg13235366", "cg13318914", "cg23945952")])
#son variables muy correlacionadas, son variables diferentes a las que utiliza el otro árbol

df2 <- data.frame(importancia = arbolTest$variable.importance, "variable" = c("cg00175838", "cg01890836", "cg04596071", "cg13235366", "cg13318914", "cg23945952"))
ggplot2::ggplot(df2) + labs(tilte="Importancia de las variables") +
  geom_segment(aes(x = variable, y = importancia, xend = variable, yend = importancia), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = variable, y = importancia, col = variable), 
             size = 4, show.legend = F) +
  coord_flip() +
  theme_bw() + ggtitle("Importancia de las variables")

#miro si las variables de ambos árboles están correlacionadas
cor(btotal[,c("cg00175838", "cg01890836", "cg04596071", "cg13235366", "cg13318914", "cg23945952","cg04903825","cg05372993", "cg22956858")])>0.8

#aciertos y fallos del árbol
pred<-predict(arbolTest,bTest, type='class')
t2<- table(bTest$Group, pred)
prop.table(t2,1) #probabilidad de fallar cada tipo
t(t2) #no comete ningún error

#evaluamos en la muestra b
predb<-predict(arbolTest,b, type='class')
tb<- table(b$Group, predb)
t(tb) #5 fallos, clasifica 1 maligno como benigno y 4 BN como M
#el error 5/39=0.13

#estudio las otras variables igual de importantes
variablesT <- names(arbolTest$variable.importance)
evaluar_variables(variablesT, bTest, b)
#la mejor es cg00175838


#lo hago con la muestra total
arboltot<-rpart(Group~.,btotal,control=rpart.control(xval=67, cp=0))
print(arboltot)
printcp(arboltot)
rpart.plot(arboltot, main="", extra=2, faclen=0, type=1, yes.text="si",no.text="no")

pred<-predict(arboltot,btotal, type='class')
t<- table(btotal$Group, pred)
t(t)

arboltot$variable.importance 
cor(bTest[,c("cg04903825","cg13027371", "cg14779065", "cg16515974", "cg17548261", "cg24227782" )])
df2 <- data.frame(importancia = arboltot$variable.importance, "variable" = c("cg04903825","cg13027371", "cg14779065", "cg16515974", "cg17548261", "cg24227782"))
ggplot2::ggplot(df2) + labs(tilte="Importancia de las variables") +
  geom_segment(aes(x = variable, y = importancia, xend = variable, yend = importancia), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = variable, y = importancia, col = variable), 
             size = 4, show.legend = F) +
  coord_flip() +
  theme_bw() + ggtitle("Importancia de las variables")

variablesTotal <- names(arboltot$variable.importance)
evaluar_variables(variablesTotal, btotal, btotal)
#la mejor es cg04903825

#-----------------------------------------------------------------------------
#fuerzo un árbol más profundo con b
#para ello incluyo el parámetro minsplit=1
start_time <- Sys.time()
arbolprof<-rpart(Group~.,b,control=rpart.control(xval=39, cp=0, minsplit=1))
end_time <- Sys.time()
tiempo_arbol <- end_time - start_time
print(tiempo_arbol)
#tarda 34.5314s

print(arbolprof)
printcp(arbolprof)
plotcp(arbolprof) #no merece la pena alargarlo, no mejora el xerror, aunque acierta el 100% de la muestra b
#xerror 0.4359*0.11765=0.05

rpart.plot(arbolprof, main = "", extra=2, faclen=0, type=1, yes.text="si",no.text="no")
#sólo utiliza cg22956858 y cg07576409, la nueva variable NO aparecia antes en el top 6 de importancia 

arbolprof$variable.importance #las más importantes cg04903825 cg05372993 cg22956858 cg15625618 cg20089319 cg21137882, se repiten
arbolprofT<-rpart(Group~.,bTest,control=rpart.control(xval=28, cp=0, minsplit=1))
print(arbolprofT)
plotcp(arbolprofT)

#--------------------------------------------------------------------------------
#############################################################################################################################
#FILTRADO

#En primer lugar, empiezo utilizando un test no parametrico de wilcoxson
#no podemos garantizar TCL por tener solo 39 individuos, asi que no podemos utilizar un t.test

pvalorestot <- sapply(1:12708,
                      function(i) wilcox.test(btotal[,i]~btotal$Group)$p.value)
names(pvalorestot) <- names(btotal[1:12708])
hist(pvalorestot)
#ajusto los contrastes multiples
pvaloresajusttot<-p.adjust(pvalorestot, method="bonferroni", n=12708)
hist(pvaloresajusttot)
mean(pvalorestot<.05)*100 #el 35.16% de las cg tienen un pvalor menor que 0.05=nivel significacion
head(sort(pvalorestot),20)
cpgtot<- names(head(sort(pvaloresajusttot),2000))
cpggrouptot<- c(cpgtot,"Group")

#compruebo que siguen las del paper
match("cg17154646", names(sort(pvaloresajusttot)))
match("cg21915100", names(sort(pvaloresajusttot)))
match("cg10189462", names(sort(pvaloresajusttot))) 
match("cg26977176", names(sort(pvaloresajusttot))) 
match("cg19883872", names(sort(pvaloresajusttot))) 
#si estan

b2ktest<- bTest[,cpggrouptot]
b2k<- b[,cpggrouptot]
length(intersect(colnames(b2k),colnames(b2ktest)))
b2ktotal<- rbind(b2k, b2ktest)
#b2kTOT[, "cg00175838"]


#construyo árboles, mismo procedimiento que antes

#empiezo por la muestra bValues
arbol2k<-rpart(Group~.,b2k,control=rpart.control(xval=39, cp=0, loss=matrix(c(0,1,10,0),byrow=TRUE, nrow=2)))

#resultados: 
print(arbol2k) #cg04903825
printcp(arbol2k)
#xerror 0.15

rpart.plot(arbol2k, main = "Árbol de clasificación con la muestra bValues", extra=2, faclen=0, type=1, yes.text="si",no.text="no")

#aciertos y fallos del árbol
pred<-predict(arbol2k,b, type='class')
t<- table(b$Group, pred)
prop.table(t,1) #probabilidad de fallar cada tipo
t(t) #clasificamos 1 BN como M

#evaluamos en la muestra test
predTest<-predict(arbol2k,bTest, type='class')
tTest<- table(bTest$Group, predTest)
t(tTest) #clasifica 1M como BN
#precision = 0.96

arbol2k$variable.importance #hay tres cpg igual de importantes cg04903825 cg05372993 cg22956858
cor(b[,c("cg04903825","cg05372993", "cg22956858")])
#son variables muy correlacionadas, les siguen otras tres cpg igua de importantes cg15625618 cg20089319 cg21137882 

df <- data.frame(importancia = arbol2k$variable.importance, "variable" = c("cg04903825", "cg05372993", "cg22956858","cg00175838" ,"cg24227782", "cg26620021"))
ggplot2::ggplot(df) + labs(tilte="Importancia de las variables") +
  geom_segment(aes(x = variable, y = importancia, xend = variable, yend = importancia), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = variable, y = importancia, col = variable), 
             size = 4, show.legend = F) +
  coord_flip() +
  theme_bw() + ggtitle("Importancia de las variables")

variables2k <- names(arbol2k$variable.importance)
evaluar_variables(variables2k, b2k, b2ktest)
#la mejor es cg04903825


#estudio ahora en el caso de construirlo con la muestra test
arbol2kTest<-rpart(Group~.,b2ktest,control=rpart.control(xval=28, cp=0, loss=matrix(c(0,1,10,0),byrow=TRUE, nrow=2)))

#resultados:
print(arbol2kTest) #cg24227782
printcp(arbol2kTest)
#xerror 0

rpart.plot(arbol2kTest, main = "Árbol de clasificación con la muestra bValuesPark", extra=2, faclen=0, type=1, yes.text="si",no.text="no")

#aciertos y fallos del árbol
pred<-predict(arbol2kTest,b2ktest, type='class')
t<- table(b2ktest$Group, pred)
prop.table(t,1) #probabilidad de fallar cada tipo
t(t) #clasifica correctamente toda la muestra de entrenamiento

#evaluamos en la muestra test
predval<-predict(arbol2k,b2k, type='class')
tval<- table(b2k$Group, predval)
t(tval) #clasifica 1BN como M
#precision = 0.97

arbol2kTest$variable.importance #hay 6 cpg igual de importantes cg00175838 cg13027371 cg15797629 cg21915100 cg24227782 cg26620021 
cor(bTest[,c("cg00175838", "cg13027371", "cg15797629", "cg21915100", "cg24227782" ,"cg26620021" )])
sum(cor(bTest[,c("cg00175838", "cg13027371", "cg15797629", "cg21915100", "cg24227782" ,"cg26620021" )])>0.9)
#son variables correlacionadas

dfT <- data.frame(importancia = arbol2kTest$variable.importance, "variable" = c("cg00175838", "cg13027371", "cg15797629", "cg21915100", "cg24227782" ,"cg26620021" ))
ggplot2::ggplot(dfT) + labs(tilte="Importancia de las variables") +
  geom_segment(aes(x = variable, y = importancia, xend = variable, yend = importancia), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = variable, y = importancia, col = variable), 
             size = 4, show.legend = F) +
  coord_flip() +
  theme_bw() + ggtitle("Importancia de las variables")

variables2kTEST <- names(arbol2kTest$variable.importance)
evaluar_variables(variables2kTEST, b2ktest, b2k)
#la mejor es cg24227782


#hago un arbol con el total
arbol2ktot<-rpart(Group~.,b2ktotal,control=rpart.control(xval=67, cp=0))

#resultados:
print(arbol2ktot) #cg24227782
printcp(arbol2ktot)
#xerror 0.13

rpart.plot(arbol2ktot, main = "Árbol de clasificación con la muestra total", extra=2, faclen=0, type=1, yes.text="si",no.text="no")

#aciertos y fallos del árbol
pred<-predict(arbol2ktot,b2ktotal, type='class')
t<- table(b2ktotal$Group, pred)
prop.table(t,1) #probabilidad de fallar cada tipo
t(t) #clasifica 2 BN como M

arbol2ktot$variable.importance #hay 4 cpg igual de importantes cg24227782 cg13027371 cg14779065 cg16515974 
cor(b2ktotal[,c("cg24227782" ,"cg13027371", "cg14779065", "cg16515974" )])

df <- data.frame(importancia = arbol2ktot$variable.importance, "variable" = c("cg24227782" ,"cg13027371", "cg14779065", "cg16515974","cg21915100" ,"cg26620021" ))
ggplot2::ggplot(df) + labs(tilte="Importancia de las variables") +
  geom_segment(aes(x = variable, y = importancia, xend = variable, yend = importancia), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = variable, y = importancia, col = variable), 
             size = 4, show.legend = F) +
  coord_flip() +
  theme_bw() + ggtitle("Importancia de las variables")

variables2ktotal<- names(arbol2ktot$variable.importance)
evaluar_variables(variables2ktotal, b2ktotal, b2ktotal)
#todas igual menos cg21915100


#################################################################################################################

#calculo la tipificacion de las variables cpg

tipificacion <- c()
cg_columns <- grep("^cg", names(b2ktotal), value = TRUE)
for (col in cg_columns) {
  means <- by(b2ktotal[[col]], b2ktotal$Group, mean)
  sds <- by(b2ktotal[[col]], b2ktotal$Group, sd)
  resultado <- (means[1] - means[2]) / mean(c(sds[1], sds[2]))
  tipificacion[col] <- resultado
}
names(tipificacion) < cg_columns

tail(sort(abs(tipificacion)),20)
match("cg24227782", names(sort(abs(tipificacion))))
match("cg04903825", names(sort(abs(tipificacion))))
match("cg26620021", names(sort(abs(tipificacion))))
match("cg21915100", names(sort(abs(tipificacion))))
#las variables relevantes tienen buena tipificacion


#represento los datos en función de las dos variables mas discriminantes, segun los árboles

lista <- which(b2ktotal$Group == "M")
plot(b2ktotal$cg04903825, b2ktotal$cg24227782, 
     xlab = "CpG04903825", ylab = "CpG24227782", 
     main = "Relación entre CpG04903825 y CpG24227782",
     pch = 16, col = "blue", cex = 0.6)
points(b2ktotal$cg04903825[lista], b2ktotal$cg24227782[lista], 
       col = "red", pch = 16, cex = 0.8)
grid()
legend("topleft", legend = c("Benigno", "Maligno"),
       col = c("blue", "red"), pch = 16, bty = "n")

#pruebo con otras variables

plot(b2ktotal$cg26620021, b2ktotal$cg24227782, 
     xlab = "CpG26620021", ylab = "CpG24227782", 
     main = "Relación entre CpG26620021 y CpG24227782",
     pch = 16, col = "blue", cex = 0.6)
points(b2ktotal$cg26620021[lista], b2ktotal$cg24227782[lista], 
       col = "red", pch = 16, cex = 0.8)
grid()
legend("topleft", legend = c("Benigno", "Maligno"),
       col = c("blue", "red"), pch = 16, bty = "n")

plot(b2ktotal$cg26620021, b2ktotal$cg04903825, 
     xlab = "CpG26620021", ylab = "CpG04903825", 
     main = "Relación entre CpG26620021 y CpG04903825",
     pch = 16, col = "blue", cex = 0.6)
points(b2ktotal$cg26620021[lista], b2ktotal$cg04903825[lista], 
       col = "red", pch = 16, cex = 0.8)
grid()
legend("topleft", legend = c("Benigno", "Maligno"),
       col = c("blue", "red"), pch = 16, bty = "n")



#DEFINO LA NUEVA VARIABLE SUMA T
b2ktotal$T <- (b2ktotal$cg04903825+b2ktotal$cg24227782+b2ktotal$cg26620021)
b2k$T <- b2k$cg04903825+b2k$cg24227782+b2k$cg26620021
b2ktest$T <- b2ktest$cg04903825+b2ktest$cg24227782+b2ktest$cg26620021

ggplot(b2ktotal, aes(x = factor(Group), y =T , fill = Group)) + 
  geom_boxplot()+ labs(x = "Tipo Tumor", y = "T",fill='Tipo', title = "Variable T")


plot(b2ktotal$T, b2ktotal$cg04903825, 
     xlab = "CpG26620021", ylab = "CpG04903825", 
     main = "Relación entre CpG26620021 y CpG04903825",
     pch = 16, col = "blue", cex = 0.6)
points(b2ktotal$T[lista], b2ktotal$cg04903825[lista], 
       col = "red", pch = 16, cex = 0.8)
grid()
legend("topleft", legend = c("Benigno", "Maligno"),
       col = c("blue", "red"), pch = 16, bty = "n")


plot(b2ktotal$T, b2ktotal$cg24227782, 
     xlab = "T", ylab = "CpG24227782", 
     main = "Relación entre T y CpG24227782",
     pch = 16, col = "blue", cex = 0.6)
points(b2ktotal$T[lista], b2ktotal$cg24227782[lista], 
       col = "red", pch = 16, cex = 0.8)
grid()
legend("topleft", legend = c("Benigno", "Maligno"),
       col = c("blue", "red"), pch = 16, bty = "n")

#construyo modelos incluyendo la variable T
arbolSUMA<-rpart(Group~.,b2ktotal, control= rpart.control(cp=0, xval=67, loss=matrix(c(0,1,10,0),byrow=TRUE, nrow=2)))

#resultados: 
print(arbolSUMA) #T
printcp(arbolSUMA)
#xerror 0.12

rpart.plot(arbolSUMA, main = "Árbol de clasificación con la muestra TOTAL", extra=2, faclen=0, type=1, yes.text="si",no.text="no")
rpart.plot(arbolSUMA, extra=2, faclen=0, type=1, yes.text="si",no.text="no")

#aciertos y fallos del árbol
pred<-predict(arbolSUMA,b2ktotal, type='class')
t<- table(b2ktotal$Group, pred)
prop.table(t,1) #probabilidad de fallar cada tipo
t(t) #clasificamos 1 BN como M

arbolSUMA$variable.importance #La variable T es la mas importante, seguida de cg04903825 cg13027371 cg14779065 cg24227782 cg26620021
df <- data.frame(importancia = arbolSUMA$variable.importance, "variable" = c("T", "cg04903825" ,"cg13027371" ,"cg14779065" ,"cg24227782" ,"cg26620021"))
ggplot2::ggplot(df) + labs(tilte="Importancia de las variables") +
  geom_segment(aes(x = variable, y = importancia, xend = variable, yend = importancia), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = variable, y = importancia, col = variable), 
             size = 4, show.legend = F) +
  coord_flip() +
  theme_bw() + ggtitle("Importancia de las variables")

#ahora construyo un árbol con cada una de las muestras individualmente
arbolSUMAb<-rpart(Group~.,b2k, control= rpart.control(cp=0, xval=39, loss=matrix(c(0,1,10,0),byrow=TRUE, nrow=2)) )

#resultados: 
print(arbolSUMAb) #cg04903825
printcp(arbolSUMAb)
#xerror 0.18

rpart.plot(arbolSUMAb, main = "Árbol de clasificación con la muestra bValues", extra=2, faclen=0, type=1, yes.text="si",no.text="no")

#aciertos y fallos del árbol
pred<-predict(arbolSUMAb,b2k, type='class')
t<- table(b2k$Group, pred)
prop.table(t,1) #probabilidad de fallar cada tipo
t(t) #clasificamos 1 BN como M

#evaluamos en la muestra test
predTest<-predict(arbolSUMAb,b2ktest, type='class')
tTest<- table(b2ktest$Group, predTest)
t(tTest) #clasifica 1M como BN

arbolSUMAb$variable.importance #hay tres cpg igual de importantes cg04903825 cg05372993 cg22956858 T
cor(b2k[,c("cg04903825" ,"cg05372993" ,"cg22956858", "T")])

df <- data.frame(importancia = arbolSUMAb$variable.importance, "variable" = c("cg04903825" ,"cg05372993" ,"cg22956858", "T",  "cg24227782" ,"cg26620021"))
ggplot2::ggplot(df) + labs(tilte="Importancia de las variables") +
  geom_segment(aes(x = variable, y = importancia, xend = variable, yend = importancia), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = variable, y = importancia, col = variable), 
             size = 4, show.legend = F) +
  coord_flip() +
  theme_bw() + ggtitle("Importancia de las variables")

variables2kT <- names(arbolSUMAb$variable.importance)
evaluar_variables(variables2kT, b2k, b2ktest)
#la mejor es T,cg24227782, cg26620021


#AHORA CON BVALUESPARK
arbolSUMAbp<-rpart(Group~.,b2ktest, control= rpart.control(cp=0, xval=28, loss=matrix(c(0,1,10,0),byrow=TRUE, nrow=2)) )

#resultados: 
print(arbolSUMAbp) #cg24227782
printcp(arbolSUMAbp)
#xerror 0

rpart.plot(arbolSUMAbp, main = "Árbol de clasificación con la muestra bValuesPark", extra=2, faclen=0, type=1, yes.text="si",no.text="no")

#aciertos y fallos del árbol
predTest<-predict(arbolSUMAbp,b2ktest, type='class')
tTest<- table(b2ktest$Group, predTest)
t(tTest) #0 fallos

#evaluamos en la muestra test
pred<-predict(arbolSUMAbp,b2k, type='class')
t<- table(b2k$Group, pred)
prop.table(t,1) #probabilidad de fallar cada tipo
t(t) #clasificamos 2 BN como M

arbolSUMAbp$variable.importance #hay 6 cpg igual de importantes cg00175838 cg13027371 cg15797629 cg21915100 cg24227782 cg26620021

df <- data.frame(importancia = arbolSUMAbp$variable.importance, "variable" = c("cg00175838" ,"cg13027371", "cg15797629", "cg21915100", "cg24227782", "cg26620021"))
ggplot2::ggplot(df) + labs(tilte="Importancia de las variables") +
  geom_segment(aes(x = variable, y = importancia, xend = variable, yend = importancia), 
               size = 1.5, alpha = 0.7) +
  geom_point(aes(x = variable, y = importancia, col = variable), 
             size = 4, show.legend = F) +
  coord_flip() +
  theme_bw() + ggtitle("Importancia de las variables")

variables2ktest <- names(arbolSUMAbp$variable.importance)
evaluar_variables(variables2ktest, b2ktest, b2k)
#la mejor es cg24227782 

#parece que con T, cg24227782 y cg04903825, sería suficiente para clasificar los nodulos

arbolSUMAb_def<-rpart(Group~.,b2k[, c("cg04903825", "cg24227782", "T" , "Group")], control= rpart.control(cp=0, xval=39, loss=matrix(c(0,1,10,0),byrow=TRUE, nrow=2)) )
print(arbolSUMAb_def)
printcp(arbolSUMAb_def)
arbolSUMAb_def$variable.importance
t(table(b2k$Group, predict(arbolSUMAb_def,b2k, type='class')))
t(table(b2ktest$Group, predict(arbolSUMAb_def,b2ktest, type='class')))
rpart.plot(arbolSUMAb_def, main = "SUMA", extra=2, faclen=0, type=1, yes.text="si",no.text="no")


arbolSUMAbp_def<-rpart(Group~.,b2ktest[, c( "cg24227782", "T" , "cg04903825","Group")], control= rpart.control(cp=0, xval=28, loss=matrix(c(0,1,10,0),byrow=TRUE, nrow=2)) )
print(arbolSUMAbp_def)
printcp(arbolSUMAbp_def)
arbolSUMAbp_def$variable.importance
t(table(b2k$Group, predict(arbolSUMAbp_def,b2k, type='class')))
t(table(b2ktest$Group, predict(arbolSUMAbp_def,b2ktest, type='class')))
rpart.plot(arbolSUMAbp_def, main = "SUMA", extra=2, faclen=0, type=1, yes.text="si",no.text="no")

#se obtienen los mismos resultados


ggplot(b2ktotal, aes(x = factor(Group), y =cg24227782 , fill = Group)) + 
  geom_boxplot()+ labs(x = "Tipo Tumor", y = "cg24227782",fill='Tipo', title = "CpG 24227782")
ggplot(b2ktotal, aes(x = factor(Group), y =T , fill = Group)) + 
  geom_boxplot()+ labs(x = "Tipo Tumor", y = "Variable T",fill='Tipo', title = "Variable T")



