"gsMMD2.v"<-
function(obj.eSet, 
         memSubjects, 
         memIni,
         maxFlag=TRUE, 
         thrshPostProb=0.50, 
         geneNames=NULL, 
         alpha=0.05, 
         transformFlag=FALSE, 
         transformMethod="boxcox", 
         scaleFlag=TRUE, 
         criterion=c("cor", "skewness", "kurtosis"),
         minL=-10, 
         maxL=10, 
         stepL=0.1, 
         eps=1.0e-3, 
         ITMAX=100, 
         plotFlag=FALSE,
         quiet=TRUE)
{
  # get expression level matrix
  X<-exprs(obj.eSet)

  res<-gsMMD2.default.v(X, 
         memSubjects, 
         memIni,
         maxFlag, 
         thrshPostProb, 
         geneNames, 
         alpha, 
         transformFlag,
         transformMethod,
         scaleFlag,
         criterion,
         minL,
         maxL,
         stepL,
         eps,
         ITMAX,
         plotFlag,
         quiet)

  invisible(res)

}


"gsMMD2.default.v"<-
function(X, 
         memSubjects, 
         memIni,
         maxFlag=TRUE, 
         thrshPostProb=0.50, 
         geneNames=NULL, 
         alpha=0.05, 
         transformFlag=FALSE, 
         transformMethod="boxcox", 
         scaleFlag=TRUE, 
         criterion=c("cor", "skewness", "kurtosis"),
         minL=-10, 
         maxL=10, 
         stepL=0.1, 
         eps=1.0e-3, 
         ITMAX=100, 
         plotFlag=FALSE,
         quiet=TRUE)
{

  transformMethod<-match.arg(transformMethod, choices=c("boxcox", "log2", "log10", "log", "none"))
  criterion<-match.arg(criterion, c("cor", "skewness", "kurtosis"))

  X<-as.matrix(X)
  nGenes<-nrow(X)
  nSubjects<-ncol(X)
  nCases<-sum(memSubjects==1, na.rm=TRUE)
  nControls<-sum(memSubjects==0, na.rm=TRUE)

  if(sum(is.null(geneNames), na.rm=TRUE))
  { geneNames<-paste("gene", 1:nGenes, sep="") }

  cat("Programming is running. Please be patient...\n")
  lambda<-NA
  if(transformFlag)
  { 
    if(transformMethod!="none")
    {
      vec<-as.numeric(X)
      min.vec<-min(vec, na.rm=TRUE)
      if(min.vec<0)
      {
        cat("****** Begin Warning ******** \n")
        cat("Warning: Data contains non-positive values! To continue ",
           transformMethod, " transformation,\n")
        cat("We first perform the following transformation:\n")
        cat("x<-x+abs(min(x, na.rm=TRUE))+1\n")
        cat("****** End Warning ******** \n")
    
        X<-X+abs(min.vec)+1
      }
    }
    tmp<-transFunc.v(X, transformMethod, criterion, 
                   minL, maxL, stepL, eps, plotFlag, ITMAX=0) 
    if(transformMethod=="boxcox")
    { X<-tmp$dat 
      lambda<-tmp$lambda.avg
    }
    else {
     X<-tmp
    }
    if(!quiet)
    { cat(paste("Data transformation (", transformMethod, ") performed\n")) }
  }

  if(scaleFlag)
  {
    if(!quiet)
    { cat("Gene profiles are scaled so that they have mean zero and variance one!\n") }
    X<-t(apply(X, 1, scale, center=TRUE, scale=TRUE))

    # to avoid linear dependence of tissue samples after scaling
    # gene profiles, we delete a tissue sample.
    # We arbitrarily select the tissue sample, which has the biggest label number, 
    # from the tissue sample group that has larger size than the other 
    # tissue sample group. For example, if there are 6 cancer tissue samples 
    # and 10 normal tissue samples, we delete the 10-th normal tissue sample after scaling.

    if(nCases>nControls)
    { 
      pos<-which(memSubjects==1)
      pos2<-pos[nCases]
      X<-X[,-pos2]
      memSubjects<-memSubjects[-pos2]
    } else {
      pos<-which(memSubjects==0)
      pos2<-pos[nControls]
      X<-X[,-pos2]
      memSubjects<-memSubjects[-pos2]
    }
    nCases<-sum(memSubjects==1, na.rm=TRUE)
    nControls<-sum(memSubjects==0, na.rm=TRUE)
    nSubjects<-nCases+nControls
  }

  #cat("Programming is running. Please be patient...\n")
  # records initial parameter estimates
  tmpIni<-getPara.v(X, memSubjects, memIni, eps)
  paraIniRP<-tmpIni$para

  # records initial log-likelihood estimates
  llkhIni<-tmpIni$llkh

  # records E(z_{ij} | x_i, Psi^{(m)})
  wiMat<-matrix(0, nrow=nGenes, ncol=3)
  # Gene Selection based on EM algorithm
  res<- paraEst.v(X, paraIniRP, memSubjects=memSubjects, 
                 maxFlag=maxFlag, thrshPostProb, geneNames=geneNames,
                 ITMAX=ITMAX, eps=eps, quiet=quiet)

  if(res$loop==0)
  {
    paraRP<-paraIniRP
    llkh<-llkhIni
    memGenes<-memIni
    wiMat<-res$memMat
  } else {
    paraRP<-res$para
    llkh<-res$llkh
    memGenes<-res$memGenes
    wiMat<-res$memMat
  }
  para<-paraRPConverter(paraRP, nCases, nControls)
  names(para)<-paraNames
  names(paraRP)<-paraNamesRP

  names(memGenes)<-geneNames
  rownames(wiMat)<-geneNames
  colnames(wiMat)<-paste("cluster", 1:3, sep="")

  memGenes2<-rep(1, nGenes)
  memGenes2[memGenes==2]<-0 # non-differentially expressed genes

  paraIni<-paraRPConverter(paraIniRP, nCases, nControls)
  names(paraIni)<-paraNames
  names(paraIniRP)<-paraNamesRP

  if(!quiet)
  {
    cat("*******************************************************\n\n")
    cat("Initial parameter estimates>>\n"); print(round(paraIni,3)); cat("\n");
    cat("Initial loglikelihood>>\n"); print(round(llkhIni,3)); cat("\n");
    cat("Final parameter estimates>>\n"); print(round(para,3)); cat("\n");
    cat("Final loglikelihood>>\n"); print(round(llkh,3)); cat("\n");
    cat("*******************************************************\n\n")
  }
  
  res<-list(dat=X, memSubjects=memSubjects, 
            memGenes=memGenes, memGenes2=memGenes2, 
            para=para, 
            llkh=llkh, wiMat=wiMat, 
            memIni=memIni, paraIni=paraIni, llkhIni=llkhIni,
            lambda=lambda)
  invisible(res) 
}

