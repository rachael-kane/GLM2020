#' Genome-wide association analysis using a generalized linear model.
#'
#' @param y A numeric matrix containing phenotype data. Dimensions are n individuals (rows) by m genetic markers (columns).
#' @param G A numeric matrix containing genotype data. Dimensions are n individuals (rows) by 1 column.
#' @param C An optional numeric matrix of dimensions containing user-specified cofactors. Dimensions are n individuals (rows) by c cofactors (columns).
#' @param PC An integer specifying the number of cofactors to retain for analysis.
#'
#' @return A numeric matrix containing p-values for each genetic marker. Dimensions are 1 row by m genetic markers (columns).

GWASbyGLM<-function(y, G, C, PC){

  n=nrow(G)
  m=ncol(G)
  mean.Y<-mean(y)
  my<-matrix(1, nrow=n, ncol=1)


  pca.obj<-prcomp(G)
  pca<-pca.obj$x

  if(missing(C)){
    C.new<-pca[,1:PC]
  }else{
    pca.c.corr.test<-corr.test(x=C[,1:ncol(C)], y=pca[,1:ncol(pca)], adjust="none")
    sig.pca.c.corr<-pca.c.corr.test$p<0.05
    filtered.pca.temp<-matrix(ncol=1,nrow=nrow(pca))
    for (i in ncol(sig.pca.c.corr)){
      if (sum(sig.pca.c.corr[,i])==0){
        filtered.pca.temp<-cbind(filtered.pca.temp, pca[,i])
      }
    }
    filtered.pca<-data.matrix(filtered.pca.temp[,2:ncol(filtered.pca.temp)])
    retained.pca<-filtered.pca[,1:PC]

    C.new<-cbind(C,filtered.pca)
  }


  P=matrix(NA,1,m)
  for (i in 1:m){
    x=G[,i]
    if(max(x)==min(x)){
      p=1}else{
        X=cbind(my,C.new,x)
        LHS=t(X)%*%X
        C=solve(LHS)
        RHS=t(X)%*%y
        b=C%*%RHS
        yb=X%*%b
        e=y-yb
        n=length(y)
        ve=sum(e^2)/(n-1)
        vt=C*ve
        t=b/sqrt(diag(vt))
        p=2*(1-pt(abs(t),n-2))
      }
    P[,i]=p[length(p)]
  }
  return(P)
}
