treinaperceptron<- function(xin, yd, eta, tol, maxepocas, par){
  #xin: entrada Nxn de dados da matriz
  #yd: tem que ser gerado para as xin (concatenado xall), metade 0 e metade 1
  #eta: peso da atualização do passo
  #tol tolerância do erro
  #maxepocas: numero máximo de epocas permitido
  #par: par==1 indica que -1 precisa ser acrescido a xin
  
  dimxin<- dim(xin)
  N<- dimxin[1]
  n<- dimxin[2]
  
  if(par==1){
    wt<- as.matrix(runif(n+1)-0,5) #inicializa o vetor de pesos w 
    xin<-cbind(1,xin) #coloca 1 no final do vetor (represeta theta)
  }
  else wt<- as.matrix(runif(n)-0,5)
  
  nepocas<- 0
  eepocas<- tol+1
  
  evec<-matrix(nrow = 1, ncol=maxepocas) #inicializa o vetor de erro 
  while((nepocas<maxepos) && (nepocas>tol)){
    
    ei2<-0
    xseq<-sample(N) # a função sample gera uma sequência aleatória para treinamento
    
    for(i in 1:N){
      irand<-xseq[i] #pega um valor do vetor aleatório
      yhati<-1.0*((xin[irand,] %*% wt))
      ei<-yd[irand]-yhati # calcula o erro
      dw<-eta*ei*xin[irand,]
      wt<-wt+dw
      ei2<-ei2+ei*ei
    }
    nepocas<-nepocas+1
    evec[nepocas]<-ei2/N
    
    eepoca<-evec[nepocas]
    
  }
  
  retlist<-list(wt, evec[1:nepocas])
  return (retlist)
  
}