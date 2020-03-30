#' Genome-wide association analysis using a general linear model.
#'
#' @param y A numeric matrix containing phenotype data. Dimensions are n rows (individuals) by 1 column.
#' @param G A numeric matrix containing genotype data. Dimensions are n rows (individuals) by m columns (genetic markers).
#' @param C A numeric matrix containing covariate data. Dimensions are n rows (individuals) by t columns (covariates).
#'  The expected input for this parameter is the $cov numeric matrix returned from the cofactor.pca.cor function included in this package.
#' @param NC An integer specifying the number of covariates to retain for analysis.
#'
#' @return A numeric matrix containing a p-value for each genetic marker. Dimensions are 1 row by m columns (genetic markers).
#'
#' @details Numeric matrices should contain only phenotype/genotype/covariate values, no accessory information like taxa ID.

GWASbyGLM<-function(y, G, C, NC){

  #define general parameters that will be used throughout
  n=nrow(G)
  m=ncol(G)

  #establishes the mean y value as 1 (corresponding to b0)
  my<-matrix(1, nrow=n, ncol=1)

  #Creates a new matrix of cofactors by subsetting C based on the number of cofactors (NC) that the user wants to retain
  C.new<-C[,1:NC]

  #Defines a new matrix, which will be filled with p-values for each of the markers
  P=matrix(NA,1,m)

  #Loop through each marker
  for (i in 1:m){
    #Isolates the column containing that marker
    x=G[,i]
    #Identifies markers that are invariant, and assigns them the maximum p-value of 1
    if(max(x)==min(x)){
      p=1}else{
        #binds the matrices my (mean y), C.new (cofactors), and x (the genetic marker of interest)
        X=cbind(my,C.new,x)
        #matrix multiplication of X by its transpose; yields left hand side
        LHS=t(X)%*%X
        #solve equation system
        C=solve(LHS)
        #right hand side
        RHS=t(X)%*%y
        #calculate b
        b=C%*%RHS
        yb=X%*%b
        e=y-yb
        n=length(y)
        ve=sum(e^2)/(n-1)
        vt=C*ve
        t=b/sqrt(diag(vt))
        #Determines p-value of genetic marker of interest
        p=2*(1-pt(abs(t),n-2))
      }
    #Assigns the p-value of this genetic marker into the corresponding cell of the matrix being populated with p-values
    P[,i]=p[length(p)]
  }
  #Outputs the matrix of p-values
  return(P)
}
