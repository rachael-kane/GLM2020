#' Test for correlations between user-specified cofactors and principal components calculated from genotype data.
#'
#' @param U A numeric matrix containing user-specified cofactors. Dimensions are n rows (individuals) by c columns (cofactors).
#' @param G A numeric matrix containing genotype data. Dimensions are n rows (individuals) by m columns (genetic markers).
#'
#' @return A numeric matrix containing user-specified cofactors and all principal components not correlated with the
#'  user-specified cofactors. Dimensions are n rows (individuals) by c columns (cofactors).
#'
#' If U is unspecified, cofactor.pca.cor will carry out principal components analysis identically to the native R function prcomp(),
#'  and cofactor.pca.cor will return principal components scores in a numeric matrix with n rows (individuals) by c columns (principal components).

cofactor.pca.cor<-function(U, G){
  pca.obj<-prcomp(G)
  pca<-pca.obj$x

  if(missing(U)){
    gwas.covariates<-pca
  }else{
    pca.c.corr.test<-corr.test(x=U[,1:ncol(C)], y=pca[,1:ncol(pca)], adjust="none")
    sig.pca.c.corr<-pca.c.corr.test$p<0.05
    filtered.pca.temp<-matrix(ncol=1,nrow=nrow(pca))
    for (i in ncol(sig.pca.c.corr)){
      if (sum(sig.pca.c.corr[,i])==0){
        filtered.pca.temp<-cbind(filtered.pca.temp, pca[,i])
      }
    }
    filtered.pca<-data.matrix(filtered.pca.temp[,2:ncol(filtered.pca.temp)])
    gwas.covariates<-cbind(U,filtered.pca)
  }
  return(gwas.covariates)
}
