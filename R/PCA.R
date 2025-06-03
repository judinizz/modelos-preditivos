rm(list = ls())
library('RnavGraphImageData')
require(caret)
library(mvtnorm)



data(faces)
faces<-t(faces)

MostraImagem <- function(x){
  rotate<-function(x) t(apply(x,2,rev))
  Img<-matrix(x, nrow=64)
  cor<-rev(gray(50:1/50))
  Image(rotate(Img),col:cor)
}

y<-NULL
for(i in 1:nrow(faces)) {
  y<-c(y, ((i-1) %/% 10)+1)
}


X<-faces

meanx<-colMeans(X)

Xs<-X-t(replicate(dim(X)[1],meanx))

check_sanidade<- colMeans(Xs)

S<-cov(Xs)

eigS<-eigen(S)
plot(eigS$values, type ='b', xlab = 'Eixo PCA', ylab = 'Autovalor')


projX<-Xs %*% eigS$vectors

plot(projX[,1],projX[,2],type='p',xlim=c(-4,4),ylim=c(-2,2),xlab='PCA1',ylab='PCA2')

par(new=TRUE)

plot(projX[(1:50),1],projX[(1:50),2],type='p',xlim=c(-4,4),ylim=c(-2,2),col='red',xlab='PCA1',ylab='PCA2')

par(new=TRUE)

plot(projX[(51:100),1],projX[(51:100),2],type='p',xlim=c(-4,4),ylim=c(-2,2),col='blue',xlab='PCA1',ylab='PCA2')


library(MASS) 
library(caret)

set.seed(123)

classe_escolhida <- 1
y_bin <- ifelse(y == classe_escolhida, 1, 0)

meanx <- colMeans(X)
Xs <- X - matrix(meanx, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)

# PCA
S <- cov(Xs)
eigS <- eigen(S)
projX <- Xs %*% eigS$vectors


var_exp <- cumsum(eigS$values) / sum(eigS$values)
k_min <- 7 #análise gráfico autovalores


X_reduced <- projX[, 1:k_min]

# Avaliação com repetição
accuracies <- c()
repeticoes <- 10
perc_treino <- 0.6

for (i in 1:repeticoes) {
  idx <- createDataPartition(y_bin, p = perc_treino, list = FALSE)
  
  
  X_train <- X_reduced[idx, ]
  y_train <- y_bin[idx]
  
  X_test <- X_reduced[-idx, ]
  y_test <- y_bin[-idx]
  
  cat("Tamanho treino classe 1:", sum(y_train == 1), "\n")
  cat("Tamanho treino classe 0:", sum(y_train == 0), "\n")
  
  
  model_1 <- cov.wt(X_train[y_train == 1, ])
  model_0 <- cov.wt(X_train[y_train == 0, ])
  
  p1 <- mean(y_train == 1)
  p0 <- 1 - p1
  
  dens1 <- dmvnorm(X_test, mean = model_1$center, sigma = model_1$cov) * p1
  dens0 <- dmvnorm(X_test, mean = model_0$center, sigma = model_0$cov) * p0
  
  
  y_pred <- ifelse(dens1 > dens0, 1, 0)
  

  acc <- mean(y_pred == y_test)
  accuracies <- c(accuracies, acc)
  cat(sprintf("Execução %d: Acurácia = %.4f\n", i, acc))
}


cat(sprintf("\nMédia da acurácia: %.4f\n", mean(accuracies)))
cat(sprintf("Desvio padrão: %.4f\n", sd(accuracies)))








