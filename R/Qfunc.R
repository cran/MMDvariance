# created on June 28, 2012
#
# sum_{i=1}^{nGenes} log10(pi.1*f1(xi)+pi.2*f2(xi)+pi.3*f3(xi))

# negative of the M-step's objective function - the part relating to pi's
# Q(Psi)=sum_{i=1}^{nGenes} wi*V(pi)+\sum_{i=1}^{nGenes) wi*U(theta)

# theta is a 17 x 1 vector
# theta = (
#   s.c1[1], delta.n1[2], mu.c1[3], r.c1[4], mu.n1[5], r.n1[6],
#   s.2[7], mu.c2[8], r.c2[9], mu.n2[10], r.n2[11],
#   s.c3[12], delta.n3[13], mu.c3[14], r.c3[15], mu.n3[16], r.n3[17])
"negQFunc.v" <-
function(theta, dat, nc, nn, n, 
  xcMat, xnMat, xcTxc, xnTxn, 
  sumw1, sumw2, sumw3, 
  sumw1xcTxc, sumw1xcT1, sumw1xcT1.sq,
  sumw1xnTxn, sumw1xnT1, sumw1xnT1.sq,
  sumw2xcTxc, sumw2xcT1, sumw2xcT1.sq,
  sumw2xnTxn, sumw2xnT1, sumw2xnT1.sq,
  sumw3xcTxc, sumw3xcT1, sumw3xcT1.sq,
  sumw3xnTxn, sumw3xnT1, sumw3xnT1.sq
  )
{

  # log(variance) of expression levels for cluster 1 for diseased subjects
  s.c1<-theta[1] 
  # mean expression level for cluster 1 for diseased subjects
  mu.c1<-theta[3] 

  # modified logit of correlation among expression levels for cluster 1 for diseased subjects
  r.c1<-theta[4] 
  rho.c1<-(exp(r.c1)-1/(nc-1))/(1+exp(r.c1))

  # mean expression level for cluster 1 for normal subjects
  mu.n1<-theta[5]
  # modified logit of correlation among expression levels for cluster 1 for normal subjects
  r.n1<-theta[6]; 
  rho.n1<-(exp(r.n1)-1/(nn-1))/(1+exp(r.n1))
  # log(variance) of expression levels for cluster 1 for normal subjects
  delta.n1<-theta[2]
  s.n1<-s.c1-exp(delta.n1)

  # log(variance) of expression levels for cluster 2 for diseased subjects
  s.2<-theta[7]; 
  # mean expression level for cluster 2 for diseased subjects
  mu.c2<-theta[8]; 
  # modified logit of correlation among expression levels for cluster 2 for diseased subjects
  r.c2<-theta[9]; 
  rho.c2<-(exp(r.c2)-1/(nc-1))/(1+exp(r.c2))

  # mean expression level for cluster 2 for normal subjects
  mu.n2<-theta[10]; 
  # modified logit of correlation among expression levels for cluster 2 for normal subjects
  r.n2<-theta[11]; 
  rho.n2<-(exp(r.n2)-1/(nn-1))/(1+exp(r.n2))

  # log(variance) of expression levels for cluster 3 for diseased subjects
  s.c3<-theta[12]; 
  # mean expression level for cluster 3 for diseased subjects
  mu.c3<-theta[14]; 
  # modified logit of correlation among expression levels for cluster 3 for diseased subjects
  r.c3<-theta[15]; 
  rho.c3<-(exp(r.c3)-1/(nc-1))/(1+exp(r.c3))

  # mean expression level for cluster 3 for normal subjects
  mu.n3<-theta[16]
  # modified logit of correlation among expression levels for cluster 3 for normal subjects
  r.n3<-theta[17]; 
  rho.n3<-(exp(r.n3)-1/(nn-1))/(1+exp(r.n3))
  # log(variance) of expression levels for cluster 3 for normal subjects
  delta.n3<-theta[13]
  s.n3<-s.c3+exp(delta.n3)

  ###########
  sumw1acTac<-sumw1xcTxc-2*mu.c1*sumw1xcT1+nc*mu.c1^2*sumw1
  sumw1acT1.sq<-sumw1xcT1.sq-2*nc*mu.c1*sumw1xcT1+nc^2*mu.c1^2*sumw1
 
  sumw1anTan<-sumw1xnTxn-2*mu.n1*sumw1xnT1+nn*mu.n1^2*sumw1
  sumw1anT1.sq<-sumw1xnT1.sq-2*nn*mu.n1*sumw1xnT1+nn^2*mu.n1^2*sumw1
 
  ######
  sumw2acTac<-sumw2xcTxc-2*mu.c2*sumw2xcT1+nc*mu.c2^2*sumw2
  sumw2acT1.sq<-sumw2xcT1.sq-2*nc*mu.c2*sumw2xcT1+nc^2*mu.c2^2*sumw2
 
  sumw2anTan<-sumw2xnTxn-2*mu.n2*sumw2xnT1+nn*mu.n2^2*sumw2
  sumw2anT1.sq<-sumw2xnT1.sq-2*nn*mu.n2*sumw2xnT1+nn^2*mu.n2^2*sumw2
 
  ######
  sumw3acTac<-sumw3xcTxc-2*mu.c3*sumw3xcT1+nc*mu.c3^2*sumw3
  sumw3acT1.sq<-sumw3xcT1.sq-2*nc*mu.c3*sumw3xcT1+nc^2*mu.c3^2*sumw3
 
  sumw3anTan<-sumw3xnTxn-2*mu.n3*sumw3xnT1+nn*mu.n3^2*sumw3
  sumw3anT1.sq<-sumw3xnT1.sq-2*nn*mu.n3*sumw3xnT1+nn^2*mu.n3^2*sumw3

  ##########
  # for f1
  ##########
  part1.1<- -nc*log(2*pi)/2-nc*s.c1/2-nc*log(nc)/2+(nc-1)*log(nc-1)/2
  tt<-exp(r.c1)
  if(tt!=Inf)
  {
    part1.1<- part1.1+nc*log(1+exp(r.c1))/2-r.c1/2
  } else {
    part1.1<- part1.1+nc*r.c1/2-r.c1/2
  }

  part1.1<- part1.1 * sumw1

  part1.2<- -(nc-1)*(exp(-s.c1)+exp(r.c1-s.c1))/(2*nc)*sumw1acTac
  part1.3<- ((nc-2)*exp(-s.c1)+exp(r.c1-s.c1)*(nc-1)-exp(-r.c1-s.c1))/(2*nc^2)*sumw1acT1.sq

  part1.4<- -nn*log(2*pi)/2-nn*s.n1/2-nn*log(nn)/2+(nn-1)*log(nn-1)/2
  tt<-exp(r.n1)
  if(tt!=Inf)
  {
    part1.4<- part1.4+nn*log(1+exp(r.n1))/2-r.n1/2
  } else {
    part1.4<- part1.4+nn*r.n1/2-r.n1/2
  }
  part1.4<- part1.4 * sumw1

  tt<-exp(r.n1-s.n1)
  if(tt!=Inf)
  { 
    part1.5<- -(nn-1)*(exp(-s.n1)+exp(r.n1-s.n1))/(2*nn)*sumw1anTan
    part1.6<- ((nn-2)*exp(-s.n1)+exp(r.n1-s.n1)*(nn-1)-exp(-r.n1-s.n1))/(2*nn^2)*sumw1anT1.sq
    part1<-part1.1+part1.2+part1.3+part1.4+part1.5+part1.6
  } else {
    part1<-part1.1+part1.2+part1.3+part1.4
  }
#
#  part1<-part1.1+part1.2+part1.3+part1.4+part1.5+part1.6
#

  ##########
  # for f2
  ##########
  part2.1<- -nc*log(2*pi)/2-nc*s.2/2-nc*log(nc)/2+(nc-1)*log(nc-1)/2
  tt<-exp(r.c2)
  if(tt!=Inf)
  {
    part2.1<- part2.1+nc*log(1+exp(r.c2))/2-r.c2/2
  } else {
    part2.1<- part2.1+nc*r.c2/2-r.c2/2
  }

  part2.1<- part2.1 * sumw2

  part2.2<- -(nc-1)*(exp(-s.2)+exp(r.c2-s.2))/(2*nc)*sumw2acTac
  part2.3<- ((nc-2)*exp(-s.2)+exp(r.c2-s.2)*(nc-1)-exp(-r.c2-s.2))/(2*nc^2)*sumw2acT1.sq

  part2.4<- -nn*log(2*pi)/2-nn*s.2/2-nn*log(nn)/2+(nn-1)*log(nn-1)/2
  tt<-exp(r.n2)
  if(tt!=Inf)
  {
    part2.4<- part2.4+nn*log(1+exp(r.n2))/2-r.n2/2
  } else {
    part2.4<- part2.4+nn*r.n2/2-r.n2/2
  }
  part2.4<- part2.4 * sumw2

  tt<-exp(r.n2-s.2)
  if(tt!=Inf)
  { 
    part2.5<- -(nn-1)*(exp(-s.2)+exp(r.n2-s.2))/(2*nn)*sumw2anTan
    part2.6<- ((nn-2)*exp(-s.2)+exp(r.n2-s.2)*(nn-1)-exp(-r.n2-s.2))/(2*nn^2)*sumw2anT1.sq
    part2<-part2.1+part2.2+part2.3+part2.4+part2.5+part2.6
  } else {
    part2<-part2.1+part2.2+part2.3+part2.4
  }
#
#  part2<-part2.1+part2.2+part2.3+part2.4+part2.5+part2.6
#

  ##########
  # for f3
  ##########
  part3.1<- -nc*log(2*pi)/2-nc*s.c3/2-nc*log(nc)/2+(nc-1)*log(nc-1)/2
  tt<-exp(r.c3)
  if(tt!=Inf)
  {
    part3.1<- part3.1+nc*log(1+exp(r.c3))/2-r.c3/2
  } else {
    part3.1<- part3.1+nc*r.c3/2-r.c3/2
  }

  part3.1<- part3.1 * sumw3

  part3.2<- -(nc-1)*(exp(-s.c3)+exp(r.c3-s.c3))/(2*nc)*sumw3acTac
  part3.3<- ((nc-2)*exp(-s.c3)+exp(r.c3-s.c3)*(nc-1)-exp(-r.c3-s.c3))/(2*nc^2)*sumw3acT1.sq

  part3.4<- -nn*log(2*pi)/2-nn*s.n3/2-nn*log(nn)/2+(nn-1)*log(nn-1)/2
  tt<-exp(r.n3)
  if(tt!=Inf)
  {
    part3.4<- part3.4+nn*log(1+exp(r.n3))/2-r.n3/2
  } else {
    part3.4<- part3.4+nn*r.n3/2-r.n3/2
  }
  part3.4<- part3.4 * sumw3

  tt<-exp(r.n3-s.n3)
  if(tt!=Inf)
  { 
    part3.5<- -(nn-1)*(exp(-s.n3)+exp(r.n3-s.n3))/(2*nn)*sumw3anTan
    part3.6<- ((nn-2)*exp(-s.n3)+exp(r.n3-s.n3)*(nn-1)-exp(-r.n3-s.n3))/(2*nn^2)*sumw3anT1.sq
    part3<-part3.1+part3.2+part3.3+part3.4+part3.5+part3.6
  } else {
    part3<-part3.1+part3.2+part3.3+part3.4
  }
#
#  part3<-part3.1+part3.2+part3.3+part3.4+part3.5+part3.6
#
  #Qfunc<-part0+part1+part2+part3
  Qfunc<-part1+part2+part3

  if(!is.na(Qfunc))
  { 
    return(-Qfunc)
  } else {
    return(1.0e+308)
  }
}


