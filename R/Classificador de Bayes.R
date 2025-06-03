library(mlbench)
library(RSNNS)

data(PimaIndiansDiabetes)
df <- PimaIndiansDiabetes
df$diabetes <- ifelse(df$diabetes == "pos", 1, -1)
x_scaled <- scale(as.matrix(df[, -9]))
y <- df$diabetes


set.seed(42)
xy <- splitForTrainingAndTest(x_scaled, y, ratio = 0.8)
xtrain <- xy$inputsTrain
ytrain <- xy$targetsTrain
xtest <- xy$inputsTest
ytest <- xy$targetsTest


kdemulti <- function(xi, xall, h){
  N <- nrow(xall)
  n <- ncol(xall)
  diffs <- sweep(xall, 2, xi)  # xall - xi
  d2 <- rowSums(diffs^2)
  dens <- sum(exp(-d2 / (2*h^2))) / (N * (sqrt(2*pi)*h)^n)
  return(dens)
}


xc1 <- xtrain[ytrain == -1, ]
xc2 <- xtrain[ytrain == 1, ]
n1 <- nrow(xc1)
n2 <- nrow(xc2)
pc1 <- n1 / (n1 + n2)
pc2 <- n2 / (n1 + n2)


h_seq <- seq(0.1, 2, by=0.1)
acc_vec <- numeric(length(h_seq))

for (i in 1:length(h_seq)) {
  h <- h_seq[i]
  
  pxc1_test <- apply(xtest, 1, function(xi) kdemulti(xi, xc1, h))
  pxc2_test <- apply(xtest, 1, function(xi) kdemulti(xi, xc2, h))
  
  y_pred <- ifelse(pxc1_test > (pc2/pc1)*pxc2_test, -1, 1)
  acc_vec[i] <- sum(y_pred == ytest) / length(ytest)
  
}

plot(h_seq, acc_vec, type="b", xlab="h", ylab="Acurácia", main="Desempenho x h")
best_h <- h_seq[which.max(acc_vec)]
cat("Melhor h:", best_h, "\n")


h <- best_h
pxc1vec <- apply(xtrain, 1, function(xi) kdemulti(xi, xc1, h))
pxc2vec <- apply(xtrain, 1, function(xi) kdemulti(xi, xc2, h))
proj <- cbind(pxc1vec, pxc2vec)

colvec <- ifelse(ytrain == -1, "red", "blue")
plot(proj[,1], proj[,2], col=colvec, pch=19, 
     xlab="p(x|c1)", ylab="p(x|c2)",
     main=paste("Espaço de Similaridades - h =", round(h, 2)))
legend("topright", legend=c("Classe -1", "Classe +1"), col=c("red", "blue"), pch=19)

plot(log(pxc1c2[,1]), log(pxc1c2[,2]), col=colvec[1+(yall+1)/2], 
     xlab="log(pxc1)", ylab="log(pxc2)")


