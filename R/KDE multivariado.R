library('RSNNS')

kdemulti <- function(xi,xall,h){
  N<-dim(xall)[1]
  n<-dim(xall)[2]
  
  xirow<-matrix(xi, ncol=n,nrow=1)
  xirep<-matrix(xirow,ncol=n,nrow=N, byrow=T)
  matdif<-(xall-xirep)*(xall-xirep)
  dximat<-rowSums(matdif)/(h*h)
  emat<-(exp(-dximat/2))
  pxi<-sum(emat)/(N*(sqrt(2*pi)*h)^n)
  
  return(pxi)
  
}

N<-30
h<-0.5 #superposição (Qual a característica do máximo de h até começar a cair a acurácia?)
m1<-c(2,2)
m2<-c(4,4)
m3<-c(2,4)
m4<-c(4,2)

g1<-matrix(rnorm(N*2,sd=0.5),nrow=N,ncol=2)+matrix(m1,nrow=N,ncol=2,byrow=T)
g2<-matrix(rnorm(N*2,sd=0.5),nrow=N,ncol=2)+matrix(m2,nrow=N,ncol=2,byrow=T)
xc1<- rbind(g1,g2)

g3<-matrix(rnorm(N*2,sd=0.5),nrow=N,ncol=2)+matrix(m3,nrow=N,ncol=2,byrow=T)
g4<-matrix(rnorm(N*2,sd=0.5),nrow=N,ncol=2)+matrix(m4,nrow=N,ncol=2,byrow=T)
xc2<-rbind(g3,g4)


xall_1<-rbind(g1,g2,g3,g4)
yall_1<-rbind(matrix(-1,nrow=(2*N),ncol=1),matrix(1,nrow=(2*N),ncol=1))

xy<-splitForTrainingandTest(xall_1,yall_1,ratio=0.1)

xall<- xy$inputsTrain
yall<-xy$targetsTrain

idxc1<-which(yall== -1)
idxc2<-which(yall==1)
xc1_tr<-xall[idxc1,]
xc2_tr<-xall[idxc2,]

nall_tst<-dim(xtst)[1]
nall<-dim(xall)[1]
n1<-dim(xc1)[1]
n2<-dim(xc2)[1]

xtst<-xy$inputsTest
ytst<-xy$targetsTest

pc1<-n1/nall
pc2<-n2/nall

pxc1<-kdemulti(xc1[1,],xc1,h)
pxc1<-kdemulti(xc2[1,],xc2,h)

seq1x2<-seq(0,6,0.1)
lseq<-length(seq1x2)

yhat_tr<-matrix(nrow=nall,ncol=2)
yhat_tst<-matrix(nrow=nall_tst,ncol=2)

MZ<-matrix(nrow=lseq,ncol=lseq)
cr<-0

for (i in 1:lseq){
  for (j in 1:lseq){
    cr<-cr+1
    x1<-seq1x2[i]
    x2<-seq1x2[j]
    x1x2<-as.matrix((cbind(x1,x2)))
    pxc1<-kdemulti(x1x2,xc1,h)
    pxc2<-kdemulti(x1x2,xc2,h)
    MZ[i,j]<-1*(pxc1>(pc2/pc1)*pxc2)
    
      }
}

contour(seq1x2,seq1x2,MZ,xlim=c(0,6),ylim=c(0,6),xlab="x1", ylab="x2",nlevels=1)
par(new=T)
plot(xc1[,1],xc1[,2],col='red', xlim=c(0,6), ylim=c(0,6), xlab="x1", ylab="x2")
par(new=T)
plot(xc2[,1],xc2[,2],col='blue', xlim=c(0,6), ylim=c(0,6), xlab="x1", ylab="x2")

persp3d(seq1x2,seq1x2,MZ,nlevels=1,xlim=c(0,6), ylim=c(0,6),xlab="x1", ylab="x2", col='red' )


#Espaço de verossimilhanças 

pxc1vec<-matrix()
pxc2vec<-matrix()

for (i in 1:nall)
  {
  pxc1vec[i]<-kdemulti(xall[i,],xc1,h)
  pxc2vec[i]<-kdemulti(xall[i,],xc2,h)
  
}

pxc1c2<-cbind(pxc1vec,pxc2vec)
colvec<-c('red','blue')

plot(pxc1c2[,1], pxc1c2[,2],col=colvec[1+(yall+1)/2], xlab="pxc1", ylab="pxc2")



