#KNN


gera_gaussianas_2classes_2D_teste <- function(N1,M1,SD1,N2,M2,SD2,seed){
  #N1 = número de amostras classe 1
  #N2 = número de amostras classe 2
  #M1 = vetor com o ponto médio da distribuição da classe 1. Ex: c(2,2) Para centrar no ponto 2,2
  #M2 = vetor com o ponto médio da distribuição da classe 2.
  #SD1 = dispersão das amostras da classe 1 (a mesma para as duas dimensões)
  #SD2 = dispersão das amostras da classe 2 (a mesma para as duas dimensões)
  #seed = semente pra se gerar o mesmo conjunto de dados
  set.seed(seed)
  
  xc1<-matrix(rnorm(N1*2),ncol=2)*SD1 + t(matrix(M1,ncol=N1,nrow=2)) #rnorm gera número aleatórios de uma distribuição normal
  xc2<-matrix(rnorm(N2*2),ncol=2)*SD2 + t(matrix(M2,ncol=N2,nrow=2)) #a matriz transposta é para ponto médio da distribuição
  
  y1 = array(1,c(N1,1))  #c gera um vetor (Tamanho do array) os elementos de c são as dimensões (N1 linhas e 1 coluna)
  y2 = array(-1,c(N2,1)) 
  
  X = rbind(xc1,xc2)  #rbind combina linhas (junta as duas matrizes)
  Y = rbind(y1,y2)
  
  
  retlist<-list(X,Y)
  
  return(retlist)  
  
}  

gera_gaussianas_2classes_2D <- function(N1, M1, SD1, N2, M2, SD2, seed) {
  set.seed(seed)
  xc1 <- matrix(rnorm(N1 * 2, mean = 0, sd = 1), ncol = 2) * SD1 + matrix(rep(M1, N1), ncol = 2, byrow = TRUE)
  
  xc2 <- matrix(rnorm(N2 * 2, mean = 0, sd = 1), ncol = 2) * SD2 + matrix(rep(M2, N2), ncol = 2, byrow = TRUE)
  
  y1 <- array(1, c(N1, 1)) 
  y2 <- array(-1, c(N2, 1)) 
  
  X <- rbind(xc1, xc2)
  Y <- rbind(y1, y2)
  
  return(list(X, Y))
}


gera_gaussianas_2classes_2D(100,c(2,2),0.941,100,c(4,4),0.941,7)

euclidean_dist <- function(x, y) sqrt(sum((x - y)^2))


knn <- function(X_train, Y_train, v, k) {
  distances <- apply(X_train, 1, function(row) euclidean_dist(row, v)) #apply aplica uma mesma função para linhas ou colunas de uma matriz, 

  neighbor_indices <- order(distances)[1:k]

  neighbor_labels <- Y_train[neighbor_indices]

  predicted_class <- ifelse(sum(neighbor_labels) >= 0, 1, -1)
  
  return(predicted_class)
}

data <- gera_gaussianas_2classes_2D(100, c(2,2), 0.941, 100, c(4,4), 0.941, 7)
X_train <- data[[1]]
Y_train <- data[[2]]
                
seqi<-seq(0.06,6,0.06) #seq(from, to, by= )
seqj<-seq(0.06,6,0.06)
M1 <- matrix(0,nrow=length(seqi),ncol=length(seqj)) 
ci<-0
for (i in seqi){
  ci<-ci+1
  cj<-0
  for(j in seqj)
  {
    cj<-cj+1
    
    M1[ci, cj] <- knn(X_train, Y_train, c(i, j), 5)
    
  }
}



library(ggplot2)

df_grid <- expand.grid(x = seqi, y = seqj)
df_grid$class <- as.vector(M1)

df_train <- as.data.frame(X_train)
df_train$label <- as.factor(Y_train)

ggplot() +
  geom_tile(data = df_grid, aes(x = x, y = y, fill = factor(class)), alpha = 0.3) +
  geom_point(data = df_train, aes(x = V1, y = V2, color = label), size = 2) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal()
