library(mlbench)
library(kernlab)
library(caret)


data(Glass)
dados <- Glass

classe_majoritaria <- names(which.max(table(dados$Type)))
dados$ClasseBinaria <- ifelse(dados$Type == classe_majoritaria, 1, -1)
dados$Type <- NULL


dados[1:9] <- scale(dados[1:9])


valores_C <- c(0.1, 1, 10)
valores_sigma <- c(0.01, 0.05, 0.1, 0.5)
grid <- expand.grid(C = valores_C, sigma = valores_sigma)


set.seed(123)
folds <- createFolds(dados$ClasseBinaria, k = 10)


resultados <- data.frame(C = numeric(), sigma = numeric(), acuracia_media = numeric(), desvio = numeric())


for (i in 1:nrow(grid)) {
  C_val <- grid$C[i]
  sigma_val <- grid$sigma[i]
  
  acuracias <- c()
  
  for (j in 1:10) {
    teste_idx <- folds[[j]]
    treino <- dados[-teste_idx, ]
    teste <- dados[teste_idx, ]
    
    modelo <- ksvm(as.matrix(treino[, 1:9]), as.factor(treino$ClasseBinaria),
                   type = "C-bsvc", kernel = "rbfdot",
                   kpar = list(sigma = sigma_val), C = C_val)
    
    pred <- predict(modelo, as.matrix(teste[, 1:9]))
    acc <- mean(pred == as.factor(teste$ClasseBinaria))
    acuracias <- c(acuracias, acc)
  }
  
  media_acc <- mean(acuracias)
  desvio_acc <- sd(acuracias)
  
  resultados <- rbind(resultados,
                      data.frame(C = C_val, sigma = sigma_val,
                                 acuracia_media = media_acc, desvio = desvio_acc))
}


resultados_ordenados <- resultados[order(-resultados$acuracia_media), ]
print(resultados_ordenados)


melhor <- resultados_ordenados[1, ]
cat(sprintf("\nMelhor configuração: C = %g, sigma = %g\n", melhor$C, melhor$sigma))
cat(sprintf("Acurácia média = %.4f, Desvio padrão = %.4f\n", melhor$acuracia_media, melhor$desvio))



X <- as.matrix(dados[, 1:9])
y <- dados$ClasseBinaria
modelo_final <- ksvm(X, y, type = "C-bsvc", kernel = "rbfdot",
                     kpar = list(sigma = melhor$sigma), C = melhor$C)


d <- as.matrix(dist(X))
h <- 1 / sqrt(2 * melhor$sigma)
K <- exp(- (d^2) / (2 * h^2))  # matriz kernel RBF

K1 <- K[, y == 1]
K_1 <- K[, y == -1]
sim_com_1 <- rowMeans(K1)
sim_com_m1 <- rowMeans(K_1)
proj <- cbind(sim_com_1, sim_com_m1)  # matriz de projeção 

sv_idx <- SVindex(modelo_final)

plot(proj[y == 1, 1], proj[y == 1, 2], col = "blue", pch = 17,
     xlab = "Similaridade com Classe +1", ylab = "Similaridade com Classe -1",
     main = "Projeção no Espaço de Similaridades")

points(proj[y == -1, 1], proj[y == -1, 2], col = "red", pch = 19)
points(proj[sv_idx, 1], proj[sv_idx, 2], pch = 1, col = "black", cex = 1.5, lwd = 2)
legend("topright",
       col = c("blue", "red", "black"), pch = c(17, 19, 1), pt.cex = c(1, 1, 1.5))


