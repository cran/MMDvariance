# modified on Sept. 28, 2011
#  (1) added 'na.rm=TRUE' to function 'sum'
#
# get wij=pi.j*f.j(xi)/[pi.1*f.1(xi)+pi.2*f.2(xi)+pi.3*f.3(xi)]
# j=1,2,3

#Psi.m<- c(pi.1, pi.2, 
#          s.c1, delta.n1, mu.c1, r.c1, mu.n1, r.n1,
#          s.2, mu.c2, r.c2, mu.n2, r.n2,
#          s.c3, delta.n3, mu.c3, r.c3, mu.n3, r.n3)
#

"wiFun.v" <-
function(Xi, Psi.m, memSubjects, eps=1.0e-6)
{
  # check if parameters are in appropriate ranges
  checkPara.v(Psi.m, eps)
  nc<-sum(memSubjects==1, na.rm=TRUE)
  nn<-sum(memSubjects==0, na.rm=TRUE)
  n<-nc+nn

  # mixture proportions
  pi.1<-Psi.m[1]; pi.2<-Psi.m[2]; 
  #pi.3<-Psi.m[3];
  pi.3<-1-pi.1-pi.2

  Psi.m2<-Psi.m[-c(1:2)]
  # log(variance) of expression levels for cluster 1 for diseased subjects
  s.c1<-Psi.m2[1] 
  sigma2.c1<-exp(s.c1)
  # mean expression level for cluster 1 for diseased subjects
  mu.c1<-Psi.m2[3] 

  # modified logit of correlation among expression levels for cluster 1 for diseased subjects
  r.c1<-Psi.m2[4] 
  rho.c1<-(exp(r.c1)-1/(nc-1))/(1+exp(r.c1))

  # mean expression level for cluster 1 for normal subjects
  mu.n1<-Psi.m2[5]
  # modified logit of correlation among expression levels for cluster 1 for normal subjects
  r.n1<-Psi.m2[6]; 
  rho.n1<-(exp(r.n1)-1/(nn-1))/(1+exp(r.n1))
  # log(variance) of expression levels for cluster 1 for normal subjects
  delta.n1<-Psi.m2[2]
  s.n1<-s.c1-exp(delta.n1)
  sigma2.n1<-exp(s.n1)

  # log(variance) of expression levels for cluster 2 for diseased subjects
  s.2<-Psi.m2[7]; 
  sigma2.2<-exp(s.2)
  # mean expression level for cluster 2 for diseased subjects
  mu.c2<-Psi.m2[8]; 
  # modified logit of correlation among expression levels for cluster 2 for diseased subjects
  r.c2<-Psi.m2[9]; 
  rho.c2<-(exp(r.c2)-1/(nc-1))/(1+exp(r.c2))

  # mean expression level for cluster 2 for normal subjects
  mu.n2<-Psi.m2[10]; 
  # modified logit of correlation among expression levels for cluster 2 for normal subjects
  r.n2<-Psi.m2[11]; 
  rho.n2<-(exp(r.n2)-1/(nn-1))/(1+exp(r.n2))

  # log(variance) of expression levels for cluster 3 for diseased subjects
  s.c3<-Psi.m2[12]; 
  sigma2.c3<-exp(s.c3)
  # mean expression level for cluster 3 for diseased subjects
  mu.c3<-Psi.m2[14]; 
  # modified logit of correlation among expression levels for cluster 3 for diseased subjects
  r.c3<-Psi.m2[15]; 
  rho.c3<-(exp(r.c3)-1/(nc-1))/(1+exp(r.c3))

  # mean expression level for cluster 3 for normal subjects
  mu.n3<-Psi.m2[16]
  # modified logit of correlation among expression levels for cluster 3 for normal subjects
  r.n3<-Psi.m2[17]; 
  rho.n3<-(exp(r.n3)-1/(nn-1))/(1+exp(r.n3))
  # log(variance) of expression levels for cluster 3 for normal subjects
  delta.n3<-Psi.m2[13]
  s.n3<-s.c3+exp(delta.n3)
  sigma2.n3<-exp(s.n3)

  ###########
  if(nc>=n){
    stop("Number of cases >= Total number of patients!\n")
  }
  if(nn>=n){
    stop("Number of controls >= Total number of patients!\n")
  }

  # expression levels of the gene for diseased subjects
  Xci<-Xi[memSubjects==1]
  # expression levels of the gene for non-diseased subjects
  Xni<-Xi[memSubjects==0]

  ###
  # density for genes in cluster 1 (over-expressed)
  ###
  XiTXi.c<-sum(Xci^2, na.rm=TRUE)
  XiT1.c<-sum(Xci, na.rm=TRUE)
  aiTai.c1<-XiTXi.c-2*mu.c1*XiT1.c+nc*mu.c1^2
  aiT12.c1<-XiT1.c^2-2*nc*mu.c1*XiT1.c+nc^2*mu.c1^2

  XiTXi.n<-sum(Xni^2, na.rm=TRUE)
  XiT1.n<-sum(Xni, na.rm=TRUE)
  aiTai.n1<-XiTXi.n-2*mu.n1*XiT1.n+nn*mu.n1^2
  aiT12.n1<-XiT1.n^2-2*nn*mu.n1*XiT1.n+nn^2*mu.n1^2

  part.c1<-( aiTai.c1  - rho.c1* aiT12.c1 / (1+(nc-1)*rho.c1) )/ (sigma2.c1 * (1-rho.c1)) 
  part.n1<-( aiTai.n1  - rho.n1* aiT12.n1 / (1+(nn-1)*rho.n1) )/ (sigma2.n1 * (1-rho.n1)) 
  delta1<- (part.c1+part.n1)

  log.detSigma1<-nc*log(sigma2.c1)+((nc-1)*log(1-rho.c1)+log(1+(nc-1)*rho.c1))+
                 nn*log(sigma2.n1)+((nn-1)*log(1-rho.n1)+log(1+(nn-1)*rho.n1))

  ###
  # density for genes in cluster 2 (non-differentially expressed)
  ###
  XiTXi.c<-sum(Xci^2, na.rm=TRUE)
  XiT1.c<-sum(Xci, na.rm=TRUE)
  aiTai.c2<-XiTXi.c-2*mu.c2*XiT1.c+nc*mu.c2^2
  aiT12.c2<-XiT1.c^2-2*nc*mu.c2*XiT1.c+nc^2*mu.c2^2

  XiTXi.n<-sum(Xni^2, na.rm=TRUE)
  XiT1.n<-sum(Xni, na.rm=TRUE)
  aiTai.n2<-XiTXi.n-2*mu.n2*XiT1.n+nn*mu.n2^2
  aiT12.n2<-XiT1.n^2-2*nn*mu.n2*XiT1.n+nn^2*mu.n2^2

  part.c2<-( aiTai.c2  - rho.c2* aiT12.c2 / (1+(nc-1)*rho.c2) )/ (sigma2.2 * (1-rho.c2)) 
  part.n2<-( aiTai.n2  - rho.n2* aiT12.n2 / (1+(nn-1)*rho.n2) )/ (sigma2.2 * (1-rho.n2)) 
  delta2<- (part.c2+part.n2)

  log.detSigma2<-nc*log(sigma2.2)+((nc-1)*log(1-rho.c2)+log(1+(nc-1)*rho.c2))+
                 nn*log(sigma2.2)+((nn-1)*log(1-rho.n2)+log(1+(nn-1)*rho.n2))

  ###
  # density for genes in cluster 3 (under-expressed)
  ###
  aiTai.c3<-XiTXi.c-2*mu.c3*XiT1.c+nc*mu.c3^2
  aiT12.c3<-XiT1.c^2-2*nc*mu.c3*XiT1.c+nc^2*mu.c3^2

  aiTai.n3<-XiTXi.n-2*mu.n3*XiT1.n+nn*mu.n3^2
  aiT12.n3<-XiT1.n^2-2*nn*mu.n3*XiT1.n+nn^2*mu.n3^2

  part.c3<-( aiTai.c3  - rho.c3* aiT12.c3 / (1+(nc-1)*rho.c3) )/ (sigma2.c3 * (1-rho.c3)) 
  part.n3<-( aiTai.n3  - rho.n3* aiT12.n3 / (1+(nn-1)*rho.n3) )/ (sigma2.n3 * (1-rho.n3)) 
  delta3<- (part.c3+part.n3)

  log.detSigma3<-nc*log(sigma2.c3)+(nc-1)*log(1-rho.c3)+log(1+(nc-1)*rho.c3)+
                 nn*log(sigma2.n3)+(nn-1)*log(1-rho.n3)+log(1+(nn-1)*rho.n3)

  ##############
  tt<-c(delta1, delta2, delta3)
  names(tt)<-c("delta1", "delta2", "delta3")
  pos<-which(tt==min(tt, na.rm=TRUE))
  pos<-pos[1]
  ################


  if(pos==1) {
    # wi1=pi.1*f1/(pi.1*f1+pi.2*f2+pi.3*f3)
    #    = pi.1/(pi.1+pi.2*f2/f1 + pi.3 * f3/f1)
    # wi2=pi.2*f2/f1/[pi.1+pi.2*f2/f1 + pi.3 * f3/f1]
    # wi3=pi.3*f3/f1/[pi.1+pi.2*f2/f1 + pi.3 * f3/f1]
    # f2/f1
    # log(f2/f1)=[log(|Sigma1|)-log(|Sigma2|)+delta1-delta2]/2
    logf2df1<-(log.detSigma1-log.detSigma2+delta1-delta2)/2
    f2df1<-exp(logf2df1)
    
    # f3/f1
    # log(f3/f1)=[log(|Sigma1|)-log(|Sigma3|)+delta1-delta3]/2
    logf3df1<-(log.detSigma1-log.detSigma3+delta1-delta3)/2
    f3df1<-exp(logf3df1)

    denom.1<-(pi.1+pi.2*f2df1+pi.3*f3df1)
  
    wi1<-pi.1/denom.1
    wi2<-pi.2*f2df1/denom.1
    wi3<-pi.3*f3df1/denom.1
    
  } else if(pos==2){

    #wi2=pi.2*f2/(pi.1*f1+pi.2*f2+pi.3*f3)
    #   = pi.2/(pi.1*f1/f2+pi.2 + pi.3 * f3/f2)
    #wi1=pi.1 * f1/f2 / [pi.1*f1/f2+pi.2 + pi.3 * f3/f2]
    #wi3=pi.3 * f3/f2 / [pi.1*f1/f2+pi.2 + pi.3 * f3/f2]
    # f1/f2
    # log(f1/f2)=(log|Sigma2|-log|Sigma1|+detla2-delta1)/2
    logf1df2<-(log.detSigma2-log.detSigma1+delta2-delta1)/2
    f1df2<-exp(logf1df2)
    
    # f3/f2
    # log(f3/f2)=(log|Sigma2|-log|Sigma3|+delta2-detla3)/2
    logf3df2<-(log.detSigma2-log.detSigma3+delta2-delta3)/2
    f3df2<-exp(logf3df2)

    denom.2<-(pi.1*f1df2+pi.2+pi.3*f3df2)
  
    wi1<-pi.1*f1df2/denom.2
    wi2<-pi.2/denom.2
    wi3<-pi.3*f3df2/denom.2
  } else {
    #wi3=pi.3*f3/(pi.1*f1+pi.2*f2+pi.3*f3)
    #   = pi.3/(pi.1*f1/f3+pi.2*f2/f3 + pi.3)
    #wi1=pi.1*f1/f3 / [pi.1*f1/f3+pi.2*f2/f3 + pi.3]
    #wi2=pi.2*f2/f3 / [pi.1*f1/f3+pi.2*f2/f3 + pi.3]
    # f1/f3
    # log(f1/f3)=(log|Sigma3|-log|Sigma1|+delta3-delta1)/2
    logf1df3<-(log.detSigma3-log.detSigma1+delta3-delta1)/2
    f1df3<-exp(logf1df3)
 
    # f2/f3
    # log(f2/f3)=(log|Sigma3|-log|Sigma2|+delta3-delta2)/2
    logf2df3<-(log.detSigma3-log.detSigma2+delta3-delta2)/2
    f2df3<-exp(logf2df3)
  
    denom.3<-(pi.1*f1df3+pi.2*f2df3+pi.3)


    wi1<-pi.1*f1df3/denom.3
    wi2<-pi.2*f2df3/denom.3
    wi3<-pi.3/denom.3
  }


  res<-c(wi1, wi2, wi3)
  names(res)<-c("wi1", "wi2", "wi3")

  return(res)
}

checkPara.v<-function(para, eps=1.0e-6)
{
  if(length(para) !=TNumParaRP)
  { stop("Number of parameters is not correct!\n") }

  # mixture proportions
  pi.1<-para[1]; pi.2<-para[2]; 
  pi.3<- 1-pi.1-pi.2

  if(pi.1>=1 || pi.1<=0)
  { cat("pi.1=", pi.1, " pi.2=", pi.2, " pi.3=",pi.3,"\n") 
    stop("pi.1 should be in (0, 1)!\n") }
  if(pi.2>=1 || pi.2<=0)
  { cat("pi.1=", pi.1, " pi.2=", pi.2, " pi.3=",pi.3,"\n") 
    stop("pi.2 should be in (0, 1)!\n") }
  if(pi.3>=1 || pi.3<=0)
  { cat("pi.1=", pi.1, " pi.2=", pi.2, " pi.3=",pi.3,"\n") 
    stop("pi.3 should be in (0, 1)!\n") }
  if(abs((pi.1+pi.2+pi.3)-1)>eps)
  { cat("pi.1=", pi.1, " pi.2=", pi.2, " pi.3=",pi.3,"\n") 
    stop("pi.1+pi.2+pi.3 should be equal to 1!\n") }

}


