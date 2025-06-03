library(mlbench)
library(mclust)
library(caret)


set.seed(123)
p <- mlbench.spirals(1000, 1.5, 0.05)
dados <- as.data.frame(p$x)
colnames(dados) <- c("X1", "X2")
dados$classe <- as.factor(p$classes)
folds <- createFolds(dados$classe, k = 10, list = TRUE, returnTrain = FALSE)


pdfnvar <- function(x, m, k, n) {
  ((1 / sqrt((2 * pi)^n * det(k)))) * exp(-0.5 * t(x - m) %*% solve(k) %*% (x - m))
}


pdf_mixture <- function(x, gmm, n) {
  dens <- 0
  for (j in 1:gmm$G) {
    m <- matrix(gmm$parameters$mean[, j], ncol=1)
    k <- gmm$parameters$variance$sigma[, , j]
    w <- gmm$parameters$pro[j]
    dens <- dens + w * pdfnvar(matrix(x, ncol=1), m, k, n)
  }
  return(dens)
}


bayes_classify <- function(x, gmm1, gmm2, priors) {
  px_c1 <- pdf_mixture(x, gmm1, 2)
  px_c2 <- pdf_mixture(x, gmm2, 2)
  
  post1 <- px_c1 * priors[1]
  post2 <- px_c2 * priors[2]
  
  if (post1 > post2) return(1) else return(2)
}


resultados <- data.frame(
  Fold = 1:10,
  G1 = integer(10),
  G2 = integer(10),
  Acuracia = numeric(10)
)


for (i in 1:10) {
  teste_idx <- folds[[i]]
  treino <- dados[-teste_idx, ]
  teste <- dados[teste_idx, ]
  
  treino1 <- treino[treino$classe == 1, 1:2]
  treino2 <- treino[treino$classe == 2, 1:2]
  
  gmm1 <- Mclust(treino1, G = 10:30)
  gmm2 <- Mclust(treino2, G = 10:30)
  
  p1 <- nrow(treino1) / nrow(treino)
  p2 <- nrow(treino2) / nrow(treino)
  
  predicoes <- apply(teste[, 1:2], 1, bayes_classify, gmm1=gmm1, gmm2=gmm2, priors=c(p1, p2))
  
  resultados$G1[i] <- gmm1$G
  resultados$G2[i] <- gmm2$G
  resultados$Acuracia[i] <- mean(predicoes == as.numeric(teste$classe))
}


print(acuracias)
media <- mean(acuracias)
desvio <- sd(acuracias)

cat(sprintf("Acurácia média: %.4f\n", media))
cat(sprintf("Desvio padrão: %.4f\n", desvio))

melhor_fold <- which.max(acuracias)
teste_idx <- folds[[melhor_fold]]
treino <- dados[-teste_idx, ]

treino1 <- treino[treino$classe == 1, 1:2]
treino2 <- treino[treino$classe == 2, 1:2]

gmm1 <- Mclust(treino1)
gmm2 <- Mclust(treino2)
p1 <- nrow(treino1) / nrow(treino)
p2 <- nrow(treino2) / nrow(treino)


x_seq <- seq(-2, 2, 0.05)
y_seq <- seq(-2, 2, 0.05)
grid <- expand.grid(X1 = x_seq, X2 = y_seq)

pred_grid <- apply(grid, 1, bayes_classify, gmm1=gmm1, gmm2=gmm2, priors=c(p1, p2))
grid$classe <- as.factor(pred_grid)


library(ggplot2)
ggplot() +
  geom_tile(data=grid, aes(x=X1, y=X2, fill=classe), alpha=0.3) +
  geom_point(data=treino, aes(x=X1, y=X2, color=classe), size=1) +
  labs(title="Superfície de Separação - Melhor Fold") +
  theme_minimal()

print(resultados)


cat(sprintf("Acurácia média: %.4f\n", mean(resultados$Acuracia)))
cat(sprintf("Desvio padrão : %.4f\n", sd(resultados$Acuracia)))

melhor_fold <- which.max(resultados$Acuracia)
cat(sprintf("Melhor fold: %d com acurácia %.4f\n", melhor_fold, resultados$Acuracia[melhor_fold]))

library(dplyr)


resumo_gaussianas <- resultados %>%
  group_by(G1, G2) %>%
  summarise(
    Media_Acuracia = mean(Acuracia),
    SD_Acuracia = sd(Acuracia),
    .groups = 'drop'
  ) %>%
  arrange(desc(Media_Acuracia))

print(resumo_gaussianas)

