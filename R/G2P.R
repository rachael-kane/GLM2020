G2P=function(G,h2,alpha,NQTN,distribution){
  n=nrow(G)
  m=ncol(G)
  #Sampling QTN
  QTN.position=sample(m,NQTN,replace=F)
  SNPQ=as.matrix(G[,QTN.position])
  QTN.position
  #QTN effects
  if(distribution=="norm")
  {addeffect=rnorm(NQTN,0,1)
  }else
  {addeffect=alpha^(1:NQTN)}
  #Simulate phenotype
  effect=SNPQ%*%addeffect
  effectvar=var(effect)
  residualvar=(effectvar-h2*effectvar)/h2
  residual=rnorm(n,0,sqrt(residualvar))
  y=effect+residual
  return(list(addeffect = addeffect, y=y, add = effect, residual = residual, QTN.position=QTN.position, SNPQ=SNPQ))
}
