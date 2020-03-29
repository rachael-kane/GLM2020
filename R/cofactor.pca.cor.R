#' Correlation between cofactors and principal components.
#'
#' @description Test for correlations between user-specified cofactors and principal components calculated from genotype data.
#'
#' @param U A numeric matrix containing user-specified cofactors. Dimensions are n rows (individuals) by t columns (cofactors).
#' @param G A numeric matrix containing genotype data. Dimensions are n rows (individuals) by m columns (genetic markers).
#'
#' @return A list of 1 or 3 objects.
#'
#' \code{cofactor.pca.cor}
#' When U is unspecified, cofactor.pca.cor will return a list of 1 object.
#' With U unspecified, function will carry out principal components analysis identically to the native R function prcomp(),
#'  and cofactor.pca.cor will return principal components scores in $cov.
#' $cov is a numeric matrix containing all principal components and individual scores.
#'  Dimensions are n rows (individuals) by t columns (principal components).
#'
#' When U is specified, cofactor.pca.cor will return a list of 3 objects.
#' $orig_pc is a numeric matrix containing all original principal components and individual scores.
#' $cov is a numeric matrix containing user-specified cofactors and all principal components not correlated with the
#'  user-specified cofactors. Dimensions are n rows (individuals) by t columns (cofactors).
#' $removed is a character matrix indicating which principal components were removed.
#'
#' The $cov matrix is intended for use as the "C" argument in the GWASbyGLM function included in this package.

cofactor.pca.cor<-function(U, G){
  #Carries out principal components analysis
  pca.obj<-prcomp(G)
  #Isolates the principal component scores matrix (rows as individuals, columns as principal components)
  pca<-pca.obj$x

  #If user-specified cofactors (U) are not specified, the function returns the principal component scores matrix
  if(missing(U)){
    gwas.covariates<-pca

    #Combines the original principal components scores, the final set of covariates (U + retained principal components)
    list_cov<-list(cov=gwas.covariates)

    #Output$cov is a covariate matrix for use as the argument "C" in the GWASbyGLM function
    return(list_cov)

    #If user-specified cofactors (U) are specified, the function tests for correlations between cofactors in U and principal components
  }else{
    #Borrows the matrix correlation function from the R package "psych"
    pca.c.corr.test<-corr.test(x=U[,1:ncol(U)], y=pca[,1:ncol(pca)], adjust="none")

    #Identifies pairs of U cofactors and principal components that are significantly correlated, with a Bonferroni correction for multiple testing
    #Columns in sig.pca.c.corr are principal components, rows are U cofactors
    #The sig.pca.c.corr matrix cells contain values of 1 and 0, indicating significant correlation or lack of correlation, respectively, between the principal component and U cofactor
    sig.pca.c.corr<-pca.c.corr.test$p<(0.05/(ncol(U)*ncol(pca)))

    #Creates empy matrix, to which the retained principal components and individual scores will be attached
    filtered.pca.temp<-matrix(ncol=1,nrow=nrow(pca))

    #Creates empty dataframe, to be filled with lines indicating which principal components are removed
    removal.report.temp<-matrix(ncol=1, nrow=1)

    #When a principal component is correlated with any of the U cofactors, it is removed
    for (i in ncol(sig.pca.c.corr)){
      #Columns in sig.pca.c.corr are principal components, rows are U cofactors
      #If a principal component is uncorrelated with all U cofactors, the sum down the column equals 0
      if ((sum(sig.pca.c.corr[,i]))==0){
        filtered.pca.temp<-cbind(filtered.pca.temp, pca[,i])
      }else{
        report<-paste("Removed principal component", i)
        removal.report.temp<-rbind(removal.report.temp, report)
      }
    }

    #Creates matrix consisting of principal components and individual scores for the uncorrelated principal components
    filtered.pca<-data.matrix(filtered.pca.temp[,2:ncol(filtered.pca.temp)])

    #Binds the matrix of U cofactors with the matrix of retained principal components
    gwas.covariates<-cbind(U,filtered.pca)

    #Creates the final report of which principal components were removed
    removal.report<-removal.report.temp[2:nrow(removal.report.temp),1]

    #Combines the original principal components scores, the final set of covariates (U + retained principal components, and the removal report)
    list_origpca_retainedcov_removed<-list(orig_pc=pca, cov=gwas.covariates, removed=removal.report)

    #Function returns this set of outputs when U is specified
    #Output$cov is a covariate matrix for use as the argument "C" in the GWASbyGLM function
    return(list_origpca_retainedcov_removed)
  }
}
