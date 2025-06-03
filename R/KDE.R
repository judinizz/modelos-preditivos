library(ggplot2)


fnormal1var <- function(x, m, r) {
  (1 / (sqrt(2 * pi) * r)) * exp(-0.5 * ((x - m) / r)^2)
}

kde <- function(x, dados, h) {
  dens <- sapply(x, function(xi) {
    mean(fnormal1var(xi, dados, h))
  })
  return(dens)
}


classifica_bayes <- function(x, kde1, kde2, p1, p2) {
  post1 <- kde1 * p1
  post2 <- kde2 * p2
  ifelse(post1 > post2, 1, 2)
}

set.seed(123)
N <- 40
u1 <- 2
u2 <- 4
s1 <- 0.6


xc1 <- rnorm(N, mean=u1, sd=s1)
xc2 <- rnorm(N, mean=u2, sd=s1)

dados <- data.frame(
  x = c(xc1, xc2),
  classe = factor(c(rep(1, N), rep(2, N)))
)

teste_idx <- sample(1:(2*N), size=0.1 * 2 * N)
teste <- dados[teste_idx, ]
treino <- dados[-teste_idx, ]
xgrid <- seq(0, 6, by=0.01)


h_vals <- c(0.1, 0.4, 1.0)

for (h in h_vals) {
  dens1 <- kde(xgrid, treino$x[treino$classe == 1], h)
  dens2 <- kde(xgrid, treino$x[treino$classe == 2], h)
  
  kde_test1 <- kde(teste$x, treino$x[treino$classe == 1], h)
  kde_test2 <- kde(teste$x, treino$x[treino$classe == 2], h)
  p1 <- mean(treino$classe == 1)
  p2 <- mean(treino$classe == 2)
  
  preds <- classifica_bayes(teste$x, kde_test1, kde_test2, p1, p2)
  acc <- mean(preds == as.numeric(teste$classe))

  plot(xgrid, dens1, type="l", col="red", ylim=c(0,1),
       main=paste("Densidade para h=", h, " - Acurácia=", round(acc, 2)),
       ylab="Densidade", xlab="x")
  lines(xgrid, dens2, col="blue")
  legend("topright", legend=c("Classe 1", "Classe 2"), col=c("red", "blue"), lty=1)
}

N_menor <- 10
xc1_menor <- rnorm(N_menor, mean=u1, sd=s1)
xc2_menor <- rnorm(N_menor, mean=u2, sd=s1)

dados_menor <- data.frame(
  x = c(xc1_menor, xc2_menor),
  classe = factor(c(rep(1, N_menor), rep(2, N_menor)))
)


dens1_menor <- kde(xgrid, dados_menor$x[dados_menor$classe == 1], 0.4)
dens2_menor <- kde(xgrid, dados_menor$x[dados_menor$classe == 2], 0.4)


plot(xgrid, dens1_menor, type="l", col="red", ylim=c(0,1),
     main="Densidade com menos amostras (h=0.4)",
     ylab="Densidade", xlab="x")
lines(xgrid, dens2_menor, col="blue")
legend("topright", legend=c("Classe 1", "Classe 2"), col=c("red", "blue"), lty=1)
