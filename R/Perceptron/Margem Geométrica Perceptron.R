
set.seed(42)

rm(list=ls())


n <- 2

N1 <- 50

sd1 <- 0.2

xc1 <- matrix(rnorm(N1 * n, sd = sd1), nrow = N1, ncol = n) + matrix(c(2, 2), nrow = N1, ncol = n, byrow = TRUE)



N2 <- 50

sd2 <- 0.2

xc2 <- matrix(rnorm(N2 * n, sd = sd2), nrow = N2, ncol = n) + matrix(c(4, 4), nrow = N2, ncol = n, byrow = TRUE)



X <- rbind(xc1, xc2)

X_aug <- cbind(1, X)  # adiciona termo de bias (x0 = 1)

Y <- c(rep(-1, N1), rep(+1, N2))



plot(xc1[,1], xc1[,2], col='blue', xlim=c(0,5), ylim=c(0,5), pch=16)

points(xc2[,1], xc2[,2], col='red', pch=16)







calc_margin_braga <- function(w, X_aug, Y) {
  
  margins <- (X_aug %*% w)
  
  imax_neg = which(margins == max(margins[margins < 0]))[1]
  
  imin_pos = which(margins == min(margins[margins > 0]))[1]
  
  return(list(imax_neg, imin_pos, margins))  # margem geométrica
  
}

  


max_iter <- 1000

w <- rep(0, n + 1)  # pesos iniciados com 0 (inclui bias)

margin_history <- numeric(max_iter)

w_history <- list()



for (epoch in 1:max_iter) {
  
  updated <- FALSE
  
  index = sample(nrow(X_aug))
  
  for (i in index) {
    
    if (Y[i] * (t(w) %*% X_aug[i,]) <= 0) {
      
      w <- w + Y[i] * X_aug[i,]
      
      updated <- TRUE
      
    }
    
  }
  
  if (!updated) break  # convergiu
  
}



margins <- calc_margin_braga(w, X_aug, Y)



margem_neg <- margins[[3]][margins[[1]]]

margem_pos <- margins[[3]][margins[[2]]]



plot(xc1[,1], xc1[,2], col='blue', xlim=c(0,5), ylim=c(0,5), pch=16)

points(xc2[,1], xc2[,2], col='red', pch=16)



# Hiperplano de decisão: w0 + w1*x + w2*y = 0 → y = -(w0 + w1*x)/w2

abline(a = -w[1]/w[3], b = -w[2]/w[3], col = 'black', lwd = 2)




num_runs <- 100
best_margin <- 0
best_w <- NULL

for (run in 1:num_runs) {
  w <- rep(0, n + 1)
  
  for (epoch in 1:max_iter) {
    updated <- FALSE
    index <- sample(nrow(X_aug))
    
    for (i in index) {
      if (Y[i] * (t(w) %*% X_aug[i,]) <= 0) {
        w <- w + Y[i] * X_aug[i, ]
        updated <- TRUE
      }
    }
    
    if (!updated) break
  }
  
  # Calcular margem geométrica (mínima distância até o hiperplano)
  margins <- calc_margin_braga(w, X_aug, Y)
  margin_neg <- margins[[3]][margins[[1]]]
  margin_pos <- margins[[3]][margins[[2]]]
  
  margin_geom <- min(abs(c(margin_neg, margin_pos)))
  
  if (margin_geom > best_margin) {
    best_margin <- margin_geom
    best_w <- w
  }
}

# Visualizar melhor vetor w
plot(xc1[,1], xc1[,2], col='blue', xlim=c(0,5), ylim=c(0,5), pch=16, main='Melhor vetor w')
points(xc2[,1], xc2[,2], col='red', pch=16)
abline(a = -best_w[1]/best_w[3], b = -best_w[2]/best_w[3], col = 'green', lwd = 2)
