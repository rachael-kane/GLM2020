#' Genome-wide association analysis using a general linear model.
#'
#' @param y A numeric matrix containing phenotype data. Dimensions are n rows (individuals) by 1 column.
#' @param G A numeric matrix containing genotype data. Dimensions are n rows (individuals) by m columns (genetic markers).
#' @param C A numeric matrix containing covariate data. Dimensions are n rows (individuals) by t columns (covariates).
#'  The expected input for this parameter is the $cov numeric matrix returned from the cofactor.pca.cor function included in this package.
#' @param NC An integer specifying the number of covariates to retain for analysis.
#'
#' @return A numeric matrix containing a p-value for each genetic marker. Dimensions are 1 row by m columns (genetic markers).

GWASbyGLM<-function(y, G, C, NC){

  n=nrow(G)
  m=ncol(G)
  my<-matrix(1, nrow=n, ncol=1)

  C.new<-C[,1:NC]

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
