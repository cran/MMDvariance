# get parameter estimate based on EM algorithm
#
const.b1<-2
const.b2<-2
const.b3<-2
const.b.vec<-c(const.b1, const.b2, const.b3)

#para<- c(pi.1, pi.2, pi.3,
#         sigma2.c1, sigma2.n1, mu.c1, rho.c1, mu.n1, rho.n1,
#         sigma2.2, mu.c2, rho.c2, mu.n2, rho.n2,
#         sigma2.c3, sigma2.n3, mu.c3, rho.c3, mu.n3, rho.n3)
#

# total  20 parameters
 
"paraEst.v" <-
function(X, paraIni, memSubjects, maxFlag=TRUE, thrshPostProb=0.50, 
         geneNames=NULL, ITMAX=100, eps=1.0e-03, quiet=TRUE)
{
  # check if parameters are in appropriate ranges
  checkPara.v(paraIni, eps)

  sumb<-sum(const.b.vec, na.rm=TRUE)
  const.b1.b2<-c(const.b1, const.b2)

  n<-ncol(X) # number of patients/subjects
  nc<-sum(memSubjects==1, na.rm=TRUE)  # number of cases
  nn<-sum(memSubjects==0, na.rm=TRUE)  # number of controls
  if(nc>=n){
    stop("Number of cases >= Total number of patients!\n")
  }
  if(nn>=n){
    stop("Number of controls >= Total number of patients!\n")
  }

  # start em algorithm
  Psi.m<-paraIni
  #cat("Psi.m>>\n"); print(Psi.m); cat("\n");
  nGenes<-nrow(X)

  ##############################

  xcMat<-X[,memSubjects==1,drop=FALSE]
  xnMat<-X[,memSubjects==0,drop=FALSE]

  xcTxc<-apply(xcMat, 1, function(x) {sum(x^2, na.rm=TRUE)})
  xcT1<-apply(xcMat, 1, function(x) {sum(x, na.rm=TRUE)})

  xnTxn<-apply(xnMat, 1, function(x) {sum(x^2, na.rm=TRUE)})
  xnT1<-apply(xnMat, 1, function(x) {sum(x, na.rm=TRUE)})

  loop<-0
  while(1)
  { 
    if(!quiet)
    {
      cat("************\n")
      cat("loop=", loop, "\n")
    }

    # E-step
    # obtain weights w_{ij}(\Psi^{(m)}), i.e. E(Z_{ij} | x_i, Psi^{(m)})
    wiMat<-t(apply(X, 1, wiFun.v, Psi.m=Psi.m, memSubjects=memSubjects, eps=eps))

    w1<-as.numeric(wiMat[,1])
    w2<-as.numeric(wiMat[,2])
    w3<-as.numeric(wiMat[,3])
    sumw1<-sum(w1, na.rm=TRUE)
    sumw2<-sum(w2, na.rm=TRUE)
    sumw3<-sum(w3, na.rm=TRUE)

    if(!quiet)
    { 
      cat("sumw1=", sumw1, " sumw2=", sumw2, " sumw3=", sumw3, "\n")    
    }

    ####
    sumw1xcTxc<-sum(w1*xcTxc, na.rm=TRUE)
    sumw1xcT1<-sum(w1*xcT1, na.rm=TRUE)
    sumw1xcT1.sq<-sum(w1*xcT1^2, na.rm=TRUE)
 
    sumw1xnTxn<-sum(w1*xnTxn, na.rm=TRUE)
    sumw1xnT1<-sum(w1*xnT1, na.rm=TRUE)
    sumw1xnT1.sq<-sum(w1*xnT1^2, na.rm=TRUE)
 
    ####
    sumw2xcTxc<-sum(w2*xcTxc, na.rm=TRUE)
    sumw2xcT1<-sum(w2*xcT1, na.rm=TRUE)
    sumw2xcT1.sq<-sum(w2*xcT1^2, na.rm=TRUE)
 
    sumw2xnTxn<-sum(w2*xnTxn, na.rm=TRUE)
    sumw2xnT1<-sum(w2*xnT1, na.rm=TRUE)
    sumw2xnT1.sq<-sum(w2*xnT1^2, na.rm=TRUE)
  
    ####
    sumw3xcTxc<-sum(w3*xcTxc, na.rm=TRUE)
    sumw3xcT1<-sum(w3*xcT1, na.rm=TRUE)
    sumw3xcT1.sq<-sum(w3*xcT1^2, na.rm=TRUE)
 
    sumw3xnTxn<-sum(w3*xnTxn, na.rm=TRUE)
    sumw3xnT1<-sum(w3*xnT1, na.rm=TRUE)
    sumw3xnT1.sq<-sum(w3*xnT1^2, na.rm=TRUE)
 
    # M-step
    # mixture proportions
    piVec.new<-(c(sumw1, sumw2)+const.b1.b2-1)/(nGenes+sumb-3)
    names(piVec.new)<-c("pi.1", "pi.2")

    thetaIni<-Psi.m[-c(1:2)]
    nTheta<-length(thetaIni)
    lower<-rep(-Inf, nTheta)
    upper<-rep(Inf, nTheta)
    res<-optim(par=thetaIni, fn=negQFunc.v, 
          nc=nc, nn=nn, n=n, xcMat=xcMat, xnMat=xnMat, 
          xcTxc=xcTxc, xnTxn=xnTxn, 
          sumw1=sumw1, sumw2=sumw2, sumw3=sumw3, 
          sumw1xcTxc=sumw1xcTxc, sumw1xcT1=sumw1xcT1, sumw1xcT1.sq=sumw1xcT1.sq,
          sumw1xnTxn=sumw1xnTxn, sumw1xnT1=sumw1xnT1, sumw1xnT1.sq=sumw1xnT1.sq,
          sumw2xcTxc=sumw2xcTxc, sumw2xcT1=sumw2xcT1, sumw2xcT1.sq=sumw2xcT1.sq,
          sumw2xnTxn=sumw2xnTxn, sumw2xnT1=sumw2xnT1, sumw2xnT1.sq=sumw2xnT1.sq,
          sumw3xcTxc=sumw3xcTxc, sumw3xcT1=sumw3xcT1, sumw3xcT1.sq=sumw3xcT1.sq,
          sumw3xnTxn=sumw3xnTxn, sumw3xnT1=sumw3xnT1, sumw3xnT1.sq=sumw3xnT1.sq,
          method = "L-BFGS-B",
          lower = lower, upper = upper)

    if(res$convergence)
    {
      cat("convergence=", res$convergence, "\n")
      cat("message=", res$message, "\n")
      cat("number of calls to fn and gr>>\n"); print(res$counts); cat("\n");
      cat("thetaIni>>\n"); print(thetaIni); cat("\n");
    }
 
    negQ<-res$value
    theta<-res$par

    # update parameters
    Psi.m.new<-c(piVec.new, theta)
    names(Psi.m.new)<-paraNamesRP

    if(!quiet)
    { 
      cat("-negQ = ", -negQ, "\n")
      cat("Psi.m.new>>\n"); print(round(Psi.m.new,3)); cat("\n") 
    }
      
    err.len<-sum(abs(Psi.m.new-Psi.m)<eps, na.rm=TRUE)
    if(err.len==length(Psi.m))
    {
      # the following statement is to 
      # make sure eta.new is an update from current wMat
      # so that d Q / d eta = 0
      Psi.m<-Psi.m.new
      break
    }

    Psi.m<-Psi.m.new
    loop<-loop+1

    if(loop>ITMAX)
    { cat("*********************\n") 
      cat("Warning! Number of looping (= ", loop, 
        ") exceeds ITMAX (=", ITMAX, ")!\n")
      cat("EM algorithm did not converge!\n")
      cat("*********************\n") 
      break
    }
  }

  
  if(!quiet)
  { cat("Total iterations for EM algorithm=", loop, "\n") }

  # gene cluster membership
  memMat<-t(apply(X, 1, wiFun.v, Psi.m=Psi.m, memSubjects=memSubjects, eps=eps))
  colnames(memMat)<-c("cluster1", "cluster2", "cluster3")
  rownames(memMat)<-geneNames

  # update gene cluster membership
  memGenes<-apply(memMat, 1, maxPosFun.v, maxFlag=maxFlag, 
                  thrshPostProb=thrshPostProb)
  sumw1<-mean(memMat[,1], na.rm=TRUE)
  sumw2<-mean(memMat[,2], na.rm=TRUE)
  sumw3<-mean(memMat[,3], na.rm=TRUE)

  #####
  sumw1xcTxc<-sum(w1*xcTxc, na.rm=TRUE)
  sumw1xcT1<-sum(w1*xcT1, na.rm=TRUE)
  sumw1xcT1.sq<-sum(w1*xcT1^2, na.rm=TRUE)
 
  sumw1xnTxn<-sum(w1*xnTxn, na.rm=TRUE)
  sumw1xnT1<-sum(w1*xnT1, na.rm=TRUE)
  sumw1xnT1.sq<-sum(w1*xnT1^2, na.rm=TRUE)
 
  #####
  sumw2xcTxc<-sum(w2*xcTxc, na.rm=TRUE)
  sumw2xcT1<-sum(w2*xcT1, na.rm=TRUE)
  sumw2xcT1.sq<-sum(w2*xcT1^2, na.rm=TRUE)
 
  sumw2xnTxn<-sum(w2*xnTxn, na.rm=TRUE)
  sumw2xnT1<-sum(w2*xnT1, na.rm=TRUE)
  sumw2xnT1.sq<-sum(w2*xnT1^2, na.rm=TRUE)
 
  #####
  sumw3xcTxc<-sum(w3*xcTxc, na.rm=TRUE)
  sumw3xcT1<-sum(w3*xcT1, na.rm=TRUE)
  sumw3xcT1.sq<-sum(w3*xcT1^2, na.rm=TRUE)
 
  sumw3xnTxn<-sum(w3*xnTxn, na.rm=TRUE)
  sumw3xnT1<-sum(w3*xnT1, na.rm=TRUE)
  sumw3xnT1.sq<-sum(w3*xnT1^2, na.rm=TRUE)
 
  ####
  llkh<- -negQFunc.v(Psi.m[-c(1:2)], X, nc, nn, n,
    xcMat, xnMat, xcTxc, xnTxn, 
    sumw1, sumw2, sumw3, 
    sumw1xcTxc, sumw1xcT1, sumw1xcT1.sq,
    sumw1xnTxn, sumw1xnT1, sumw1xnT1.sq,
    sumw2xcTxc, sumw2xcT1, sumw2xcT1.sq,
    sumw2xnTxn, sumw2xnT1, sumw2xnT1.sq,
    sumw3xcTxc, sumw3xcT1, sumw3xcT1.sq,
    sumw3xnTxn, sumw3xnT1, sumw3xnT1.sq)

  memGenes2<-rep(1, nGenes)
  memGenes2[memGenes==2]<-0

  if(sum(is.null(geneNames), na.rm=TRUE))
  {
    geneNames<-paste("gene", 1:nGenes, sep="")
  } 

  names(Psi.m)<-paraNamesRP
  names(memGenes)<-geneNames
  names(memGenes2)<-geneNames

  invisible(list(para=Psi.m, llkh=llkh, memGenes=memGenes, 
                 memGenes2=memGenes2, memMat=memMat, loop=loop))
}

# find the position of the maximum element of 'x'
maxPosFun.v<-function(x, maxFlag=TRUE, thrshPostProb=0.50)
{
  if(maxFlag)
  { pos<-which(x==max(x, na.rm=TRUE))
    pos<-pos[1]
    return(pos)
  }
    
  pos<-which(x>thrshPostProb)
  len<-length(pos)
  if(len==0)
  { 
    return(2) 
  } # non-differentially expressed

  pos<-which(x==max(x, na.rm=TRUE))

  return(pos)
}


