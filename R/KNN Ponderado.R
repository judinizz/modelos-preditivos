library(ggplot2)


pdfnvar <- function(x, m, K, n) {
  ((1 / (sqrt((2 * pi)^n * (det(K))))) * 
     exp(-0.5 * t(x - m) %*% solve(K) %*% (x - m)))
}

gera_gaussianas_2classes_2D <- function(N1, M1, SD1, N2, M2, SD2, seed) {
  set.seed(seed)
  xc1 <- matrix(rnorm(N1 * 2), ncol = 2) * SD1 + matrix(rep(M1, N1), ncol = 2, byrow = TRUE)
  xc2 <- matrix(rnorm(N2 * 2), ncol = 2) * SD2 + matrix(rep(M2, N2), ncol = 2, byrow = TRUE)
  
  y1 <- array(1, c(N1, 1))  
  y2 <- array(-1, c(N2, 1))
  
  X <- rbind(xc1, xc2)
  Y <- rbind(y1, y2)
  
  return(list(X, Y))
}

calcula_Q1_Q2 <- function(X_train, Y_train, v, k, h) {
  distances <- apply(X_train, 1, function(row) sqrt(sum((row - v)^2)))
  neighbor_indices <- order(distances)[1:k]
  
  X1 <- X_train[neighbor_indices, , drop = FALSE][Y_train[neighbor_indices] == 1, , drop = FALSE]
  X2 <- X_train[neighbor_indices, , drop = FALSE][Y_train[neighbor_indices] == -1, , drop = FALSE]
  
  K <- diag(rep(h^2, 2))
  Q1 <- sum(apply(X1, 1, function(xi) pdfnvar(xi, v, K, 2)))
  Q2 <- sum(apply(X2, 1, function(xi) pdfnvar(xi, v, K, 2)))
  
  return(c(Q1, Q2))
}

k_values <- c(3, 5, 9, 15)
h_values <- c(0.3, 0.5, 0.7, 0.9)

for (h in h_values) {
  for (k in k_values) {
    data <- gera_gaussianas_2classes_2D(100, c(2,2), 0.5, 100, c(4,4), 0.5, 7)
    X_train <- data[[1]]
    Y_train <- data[[2]]
    
    seqi <- seq(0, 6, 0.1)
    seqj <- seq(0, 6, 0.1)
    
    Q_values <- data.frame(Q1 = numeric(), Q2 = numeric(), label = factor())
    
    for (i in seqi) {
      for (j in seqj) {
        Q1_Q2 <- calcula_Q1_Q2(X_train, Y_train, c(i, j), k, h)
        Q_values <- rbind(Q_values, data.frame(Q1 = Q1_Q2[1], Q2 = Q1_Q2[2], label = factor(ifelse(Q1_Q2[1] >= Q1_Q2[2], 1, -1))))
      }
    }
    
    plot <- ggplot(Q_values, aes(x = Q1, y = Q2, color = label)) +
      geom_point() +
      geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
      scale_color_manual(values = c("red", "green")) +
      theme_minimal() +
      ggtitle(paste("Espaço Q1-Q2 - k =", k, ", h =", h))
    
    print(plot)
  }
}