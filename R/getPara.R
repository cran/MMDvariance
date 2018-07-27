# created on June 28, 2012
#  (1) given expression levels, gene memberships and subject disease status,
#      get parameter estimates and  the value of log likelihood function

"getPara.v" <-
function(X, memSubjects, memGenes, eps=1.0e-6)
{
  # expression levels for genes in cluster 1 (over-variance)
  Xupp<-X[memGenes==1,, drop=FALSE] 
  # expression levels for genes in cluster 3 (under-variance)
  Xlow<-X[memGenes==3,, drop=FALSE] 
  # expression levels for genes in cluster 2 (non-differentially variance)
  Xnon<-X[memGenes==2,, drop=FALSE]

  nCases<-sum(memSubjects==1, na.rm=TRUE) # number of cases
  nControls<-sum(memSubjects==0, na.rm=TRUE) # number of controls

  nc<-nCases
  nn<-nControls
  nSubj<-nCases+nControls
  n<-nSubj
  
  n.upp<-nrow(Xupp)
  n.low<-nrow(Xlow)
  n.non<-nrow(Xnon)

  n<-nrow(X)
  memGenes2<-rep(1, n)
  memGenes2[memGenes==2]<-0

  pi.1<-n.upp/n
  pi.3<-n.low/n
  pi.2<-1-pi.1-pi.3

  Xuppc<-Xupp[,memSubjects==1, drop=FALSE]
  Xuppn<-Xupp[,memSubjects==0, drop=FALSE]

  Xnonc<-Xnon[,memSubjects==1, drop=FALSE]
  Xnonn<-Xnon[,memSubjects==0, drop=FALSE]

  Xlowc<-Xlow[,memSubjects==1, drop=FALSE]
  Xlown<-Xlow[,memSubjects==0, drop=FALSE]

  # within gene variation
  v.uppc.w<-mean(apply(Xuppc, 1, var, na.rm=TRUE), na.rm=TRUE)
  v.uppn.w<-mean(apply(Xuppn, 1, var, na.rm=TRUE), na.rm=TRUE)

  v.nonc.w<-mean(apply(Xnonc, 1, var, na.rm=TRUE), na.rm=TRUE)
  v.nonn.w<-mean(apply(Xnonn, 1, var, na.rm=TRUE), na.rm=TRUE)

  v.lowc.w<-mean(apply(Xlowc, 1, var, na.rm=TRUE), na.rm=TRUE)
  v.lown.w<-mean(apply(Xlown, 1, var, na.rm=TRUE), na.rm=TRUE)

  # between gene variation
  v.uppc.b<-var(apply(Xuppc, 1, mean, na.rm=TRUE), na.rm=TRUE)
  v.uppn.b<-var(apply(Xuppn, 1, mean, na.rm=TRUE), na.rm=TRUE)

  v.nonc.b<-var(apply(Xnonc, 1, mean, na.rm=TRUE), na.rm=TRUE)
  v.nonn.b<-var(apply(Xnonn, 1, mean, na.rm=TRUE), na.rm=TRUE)

  v.lowc.b<-var(apply(Xlowc, 1, mean, na.rm=TRUE), na.rm=TRUE)
  v.lown.b<-var(apply(Xlown, 1, mean, na.rm=TRUE), na.rm=TRUE)

  # marginal correlation
  rho.c1<-v.uppc.b/(v.uppc.b+v.uppc.w)
  rho.n1<-v.uppn.b/(v.uppn.b+v.uppn.w)

  rho.c2<-v.nonc.b/(v.nonc.b+v.nonc.w)
  rho.n2<-v.nonn.b/(v.nonn.b+v.nonn.w)

  rho.c3<-v.lowc.b/(v.lowc.b+v.lowc.w)
  rho.n3<-v.lown.b/(v.lown.b+v.lown.w)

  r.c1<-log((1+(nCases-1)*rho.c1)/((1-rho.c1)*(nCases-1)))
  r.n1<-log((1+(nControls-1)*rho.n1)/((1-rho.n1)*(nControls-1)))

  r.c2<-log((1+(nCases-1)*rho.c2)/((1-rho.c2)*(nCases-1)))
  r.n2<-log((1+(nControls-1)*rho.n2)/((1-rho.n2)*(nControls-1)))

  r.c3<-log((1+(nCases-1)*rho.c3)/((1-rho.c3)*(nCases-1)))
  r.n3<-log((1+(nControls-1)*rho.n3)/((1-rho.n3)*(nControls-1)))

  mu.c1<-mean(apply(Xuppc, 1, mean, na.rm=TRUE), na.rm=TRUE)
  mu.n1<-mean(apply(Xuppn, 1, mean, na.rm=TRUE), na.rm=TRUE)
  sigma2.c1<-v.uppc.w
  sigma2.n1<-v.uppn.w

  mu.c2<-mean(apply(Xnonc, 1, mean, na.rm=TRUE), na.rm=TRUE)
  mu.n2<-mean(apply(Xnonn, 1, mean, na.rm=TRUE), na.rm=TRUE)
  sigma2.2<-mean(apply(Xnon, 1, var, na.rm=TRUE), na.rm=TRUE)

  mu.c3<-mean(apply(Xlowc, 1, mean, na.rm=TRUE), na.rm=TRUE)
  mu.n3<-mean(apply(Xlown, 1, mean, na.rm=TRUE), na.rm=TRUE)
  sigma2.c3<-v.lowc.w
  sigma2.n3<-v.lown.w

  if(sigma2.c1<=sigma2.n1 || sigma2.c3>=sigma2.n3 )
  { 
    paraIni = rep(NA, 19)
    names(paraIni)<-paraNamesRP
    llkh<- -Inf 
  } else { 
    # s.n1 = s.c1 - exp(delta.n1)
    delta.n1<-log(sigma2.c1-sigma2.n1)
    # s.n3 = s.c3 + exp(delta.n3)
    delta.n3<-log(sigma2.n3-sigma2.c3)

    paraIni<-c(pi.1, pi.2, 
               log(sigma2.c1), delta.n1, mu.c1, r.c1, mu.n1, r.n1, 
               log(sigma2.2), mu.c2, r.c2, mu.n2, r.n2,
               log(sigma2.c3), delta.n3, mu.c3, r.c3, mu.n3, r.n3) 

    names(paraIni)<-paraNamesRP
  
    mat<-t(apply(X, 1, wiFun.v, Psi.m=paraIni, memSubjects=memSubjects, eps=eps))

    w1<-as.numeric(mat[,1])
    w2<-as.numeric(mat[,2])
    w3<-as.numeric(mat[,3])
    xcMat<-X[,memSubjects==1,drop=FALSE]
    xnMat<-X[,memSubjects==0,drop=FALSE]
 
    #####
    xcTxc<-apply(xcMat, 1, function(x) {sum(x^2, na.rm=TRUE)})
    xcT1<-apply(xcMat, 1, function(x) {sum(x, na.rm=TRUE)})
 
    xnTxn<-apply(xnMat, 1, function(x) {sum(x^2, na.rm=TRUE)})
    xnT1<-apply(xnMat, 1, function(x) {sum(x, na.rm=TRUE)})
 
    #####
    sumw1<-mean(mat[,1], na.rm=TRUE)
    sumw2<-mean(mat[,2], na.rm=TRUE)
    sumw3<-mean(mat[,3], na.rm=TRUE)
 
    ##########
    sumw1xcTxc<-sum(w1*xcTxc, na.rm=TRUE)
    sumw1xcT1<-sum(w1*xcT1, na.rm=TRUE)
    sumw1xcT1.sq<-sum(w1*xcT1^2, na.rm=TRUE)
  
    sumw1xnTxn<-sum(w1*xnTxn, na.rm=TRUE)
    sumw1xnT1<-sum(w1*xnT1, na.rm=TRUE)
    sumw1xnT1.sq<-sum(w1*xnT1^2, na.rm=TRUE)
  
    ##########
    sumw2xcTxc<-sum(w2*xcTxc, na.rm=TRUE)
    sumw2xcT1<-sum(w2*xcT1, na.rm=TRUE)
    sumw2xcT1.sq<-sum(w2*xcT1^2, na.rm=TRUE)
  
    sumw2xnTxn<-sum(w2*xnTxn, na.rm=TRUE)
    sumw2xnT1<-sum(w2*xnT1, na.rm=TRUE)
    sumw2xnT1.sq<-sum(w2*xnT1^2, na.rm=TRUE)
 
    ##########
    sumw3xcTxc<-sum(w3*xcTxc, na.rm=TRUE)
    sumw3xcT1<-sum(w3*xcT1, na.rm=TRUE)
    sumw3xcT1.sq<-sum(w3*xcT1^2, na.rm=TRUE)
  
    sumw3xnTxn<-sum(w3*xnTxn, na.rm=TRUE)
    sumw3xnT1<-sum(w3*xnT1, na.rm=TRUE)
    sumw3xnT1.sq<-sum(w3*xnT1^2, na.rm=TRUE)
  
    llkh<- -negQFunc.v(paraIni[-c(1:2)], X, nc, nn, n,
      xcMat, xnMat, xcTxc, xnTxn, 
      sumw1, sumw2, sumw3, 
      sumw1xcTxc, sumw1xcT1, sumw1xcT1.sq,
      sumw1xnTxn, sumw1xnT1, sumw1xnT1.sq,
      sumw2xcTxc, sumw2xcT1, sumw2xcT1.sq,
      sumw2xnTxn, sumw2xnT1, sumw2xnT1.sq,
      sumw3xcTxc, sumw3xcT1, sumw3xcT1.sq,
      sumw3xnTxn, sumw3xnT1, sumw3xnT1.sq)
  }

  return(list(para=paraIni, llkh=llkh))
}

