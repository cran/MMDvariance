# created on June 28, 2012
#
# define some constant
PI<-3.1415926
TNumPara<-20
paraNames<- c("pi.1", "pi.2", "pi.3",
              "sigma2.c1", "sigma2.n1", "mu.c1", "rho.c1", "mu.n1", "rho.n1",
              "sigma2.2", "mu.c2", "rho.c2", "mu.n2", "rho.n2",
              "sigma2.c3", "sigma2.n3", "mu.c3", "rho.c3", "mu.n3", "rho.n3")

# re-parametrization
TNumParaRP<-19
paraNamesRP<- c("pi.1", "pi.2", 
              "s.c1", "delta.n1", "mu.c1", "r.c1", "mu.n1", "r.n1",
              "s.2", "mu.c2", "r.c2", "mu.n2", "r.n2",
              "s.c3", "delta.n3", "mu.c3", "r.c3", "mu.n3", "r.n3")

"gsMMD.v"<-
function(obj.eSet, 
         memSubjects, 
         maxFlag=TRUE, 
         thrshPostProb=0.50, 
         geneNames=NULL, 
         alpha=0.05, 
         iniGeneMethod= "myLeveneTest", 
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

  res<-gsMMD.default.v(X, 
         memSubjects, 
         maxFlag, 
         thrshPostProb, 
         geneNames, 
         alpha, 
         iniGeneMethod,
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


"gsMMD.default.v"<-
function(X, 
         memSubjects, 
         maxFlag=TRUE, 
         thrshPostProb=0.50, 
         geneNames=NULL, 
         alpha=0.05, 
         iniGeneMethod= "myLeveneTest", 
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

  posMethod<-match(iniGeneMethod, c("myFTest", "myLeveneTest", 
    "myBFTest", "myLevene.TM", "myAWvar", 
    "myiAWvar.BF", "myiAWvar.Levene", "myiAWvar.TM"))
  tmppos<-which(is.na(posMethod)==TRUE)
  if(length(tmppos)>0)
  { msg<-paste("The initial gene partition method(s):", iniGeneMethod[tmppos], " not available!\n")
   stop(msg)
  }

  X<-as.matrix(X)
  nGenes<-nrow(X)
  nSubjects<-ncol(X)
  nMethods<-length(iniGeneMethod)
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

  # records initial parameter estimates
  paraIniMatRP<-matrix(0, nrow=TNumParaRP, ncol=nMethods)
  # records initial gene-membership estimates
  memIniMat<-matrix(0, nrow=nGenes, ncol=nMethods)

  # records initial log-likelihood estimates
  llkhIniVec<-rep(0, nMethods)

  # records parameter estimates
  paraMatRP<-matrix(0, nrow=TNumParaRP, ncol=nMethods)
  # records gene-membership estimates
  memMat<-matrix(0, nrow=nGenes, ncol=nMethods)

  # records log-likelihood estimates
  llkhVec<-rep(0, nMethods)

  cat("Programming is running. Please be patient...\n")
  # records E(z_{ij} | x_i, Psi^{(m)})
  wiArray<-array(0, c(nGenes, 3, nMethods))
  for(i in 1:nMethods)
  {
    if(!quiet)
    { cat("******** initial parameter estimates method>>", 
        iniGeneMethod[i], " *******\n") }
    tmpIni<-getIniMemGenes.v(X, memSubjects, geneNames, 
                              iniGeneMethod[i], alpha, eps=eps)
    iniGeneMethod[i]<-tmpIni$iniGeneMethod
    memIniMat[,i]<-tmpIni$memGenes
    paraIniMatRP[,i]<-tmpIni$para
    llkhIniVec[i]<-tmpIni$llkh
    ttt<-paraIniMatRP[,i]
    names(ttt)<-paraNamesRP
    if(!quiet)
    { cat("paraIniMatRP[,i]>>\n"); print(round(ttt,3)); cat("\n"); }

    # Gene Selection based on EM algorithm
    res<- paraEst.v(X, tmpIni$para, memSubjects=memSubjects, 
                   maxFlag=maxFlag, thrshPostProb, geneNames=geneNames,
                   ITMAX=ITMAX, eps=eps, quiet=quiet)

    if(res$loop==0)
    {
      paraMatRP[,i]<-paraIniMatRP[,i]
      llkhVec[i]<-llkhIniVec[i]
      memMat[,i]<-memIniMat[,i]
      wiArray[,,i]<-res$memMat
    } else {
      paraMatRP[,i]<-res$para
      llkhVec[i]<-res$llkh
      memMat[,i]<-res$memGenes
      wiArray[,,i]<-res$memMat
    }
  }
  rownames(paraIniMatRP)<-paraNamesRP
  colnames(paraIniMatRP)<-iniGeneMethod

  rownames(memIniMat)<-geneNames
  colnames(memIniMat)<-iniGeneMethod
  names(llkhIniVec)<-iniGeneMethod
  rownames(paraMatRP)<-paraNamesRP

  colnames(paraMatRP)<-iniGeneMethod
  rownames(memMat)<-geneNames
  colnames(memMat)<-iniGeneMethod
  names(llkhVec)<-iniGeneMethod
  dimnames(wiArray)<-list(geneNames,
                          paste("cluster", 1:3, sep=""),
                          iniGeneMethod)                          


  # final results
  flagPi<-rep(0, nMethods)
  paraIniMat<-matrix(0, nrow=TNumPara, ncol=nMethods)
  rownames(paraIniMat)<-paraNames
  colnames(paraIniMat)<-iniGeneMethod

  paraMat<-matrix(0, nrow=TNumPara, ncol=nMethods)
  rownames(paraMat)<-paraNames
  colnames(paraMat)<-iniGeneMethod

  for(i in 1:nMethods)
  { 
    paraIniMat[,i]<-paraRPConverter(paraIniMatRP[,i], nCases, nControls)
    paraMat[,i]<-paraRPConverter(paraMatRP[,i], nCases, nControls)
    flagPi[i]<-sum(paraMat[2,i]>paraMat[1,i] & paraMat[2,i]>paraMat[3,i], na.rm=TRUE)
  }
  tmppos<-which(flagPi==0)
  if(length(tmppos))
  { llkhVec[tmppos]<- -Inf }
  if(!quiet)
  { cat("llkhVec>>\n"); print(llkhVec); cat("\n"); }
  pos<-which(llkhVec==max(llkhVec, na.rm=TRUE))
  len<-length(pos)
  tt<-sample(1:len, 1, replace =FALSE)
  pos<-pos[tt]

  memGenes<-as.vector(memMat[,pos])
  para<-paraMat[,pos]
  paraRP<-paraMatRP[,pos]
  llkh<-llkhVec[pos]
  wiMat<-wiArray[,,pos]

  memGenes2<-rep(1, nGenes)
  memGenes2[memGenes==2]<-0 # non-differentially expressed genes

  if(!quiet)
  { cat("*******************************************************\n\n") 
    cat("Initial parameter estimates>>\n"); print(round(paraIniMat,3)); cat("\n");
    cat("Initial loglikelihood>>\n"); print(round(llkhIniVec,3)); cat("\n");
    #tmpMat<-matrix(0, nrow=TNumPara, ncol=nMethods)
    #rownames(tmpMat)<-rownames(paraMat)
    #colnames(tmpMat)<-colnames(paraMat)
    #for(i in 1:nMethods)
    #{
    #  #tmpMat[,i]<-paraMat[[i]][,1]
    #  tmpMat[,i]<-paraMat[,1]
    #}
    #cat("Final parameter estimates based on initial estimates>>\n"); print(round(tmpMat,3)); cat("\n");
    cat("Final loglikelihood based on initial estimates>>\n"); print(round(llkhVec,3)); cat("\n");
    cat("Final parameter estimates>>\n"); print(round(para,3)); cat("\n");
    cat("Final loglikelihood>>\n"); print(round(llkh,3)); cat("\n");
    cat("initial membership freq>>\n")
    print(apply(memIniMat, 2, table, useNA="ifany"))
    cat("Final membership freq>>\n")
    print(apply(memMat, 2, table, useNA="ifany"))
    cat("*******************************************************\n\n")
  }

  res<-list(dat=X, memSubjects=memSubjects, 
            memGenes=memGenes, memGenes2=memGenes2, 
            para=para, 
            llkh=llkh, wiMat=wiMat, wiArray=wiArray,
            memIniMat=memIniMat, paraIniMat=paraIniMat, llkhIniVec=llkhIniVec,
            memMat=memMat, paraMat=paraMat, llkhVec=llkhVec, lambda=lambda)
  invisible(res) 
}

getIniMemGenes.v<-function(X, memSubjects, geneNames, iniGeneMethod="myLeveneTest",
                         alpha=0.05, eps=1e-6)
{
  iniGeneMethod<-match.arg(iniGeneMethod,  
			   choices=c("myFTest", "myLeveneTest",
    "myLevene.TM", "myBFTest", "myAWvar", 
    "myiAWvar.BF", "myiAWvar.Levene", "myiAWvar.TM"))

  tmp<-iniMemGenesTestFunc.v(X, memSubjects=memSubjects, testFun=iniGeneMethod, 
                        geneNames=geneNames, alpha = alpha, eps=eps)
   
  res<-list(para=tmp$para, llkh=tmp$llkh, memGenes=tmp$memGenes, 
            memGenes2=tmp$memGenes2, iniGeneMethod=iniGeneMethod)
  return(res)
}

paraRPConverter<-function(paraRP, nCases, nControls)
{
  nc<-nCases
  nn<-nControls
  n<-nc+nn

  # mixture proportions
  pi.1<-paraRP[1]; pi.2<-paraRP[2]; 
  pi.3<-1-pi.1-pi.2

  paraRP2<-paraRP[-c(1:2)]
  #################################
  # for f1
  #################################
  # variance of expression levels for cluster 1 for diseased subjects
  s.c1<-paraRP2[1]; 
  sigma2.c1<-exp(s.c1)

  # mean expression level for cluster 1 for diseased subjects
  mu.c1<-paraRP2[3]; 

  # modified logit of correlation among expression levels for cluster 1 for diseased subjects
  r.c1<-paraRP2[4]; 
  rho.c1<-(exp(r.c1)-1/(nc-1))/(1+exp(r.c1))

  # mean expression level for cluster 1 for normal subjects
  mu.n1<-paraRP2[5]

  # modified logit of correlation among expression levels for cluster 1 for normal subjects
  r.n1<-paraRP2[6]; 
  rho.n1<-(exp(r.n1)-1/(nn-1))/(1+exp(r.n1))

  # variance of expression levels for cluster 1 for normal subjects
  delta.n1<-paraRP2[2]; 
  s.n1<-s.c1-exp(delta.n1)
  sigma2.n1<-exp(s.n1)

  #################################
  # for f2
  #################################
  # variance of expression levels for cluster 2 for diseased subjects
  s.2<-paraRP2[7]; 
  sigma2.2<-exp(s.2)

  # mean expression level for cluster 2 for diseased subjects
  mu.c2<-paraRP2[8]; 

  # modified logit of correlation among expression levels for cluster 2 for diseased subjects
  r.c2<-paraRP2[9]; 
  rho.c2<-(exp(r.c2)-1/(nc-1))/(1+exp(r.c2))

  # mean expression level for cluster 2 for normal subjects
  mu.n2<-paraRP2[10]

  # modified logit of correlation among expression levels for cluster 2 for normal subjects
  r.n2<-paraRP2[11]; 
  rho.n2<-(exp(r.n2)-1/(nn-1))/(1+exp(r.n2))

  #################################
  # for f3
  #################################
  # variance of expression levels for cluster 3 for diseased subjects
  s.c3<-paraRP2[12]; 
  sigma2.c3<-exp(s.c3)

  # mean expression level for cluster 3 for diseased subjects
  mu.c3<-paraRP2[14]; 

  # modified logit of correlation among expression levels for cluster 3 for diseased subjects
  r.c3<-paraRP2[15]; 
  rho.c3<-(exp(r.c3)-1/(nc-1))/(1+exp(r.c3))

  # mean expression level for cluster 3 for normal subjects
  mu.n3<-paraRP2[16]

  # modified logit of correlation among expression levels for cluster 3 for normal subjects
  r.n3<-paraRP2[17]; 
  rho.n3<-(exp(r.n3)-1/(nn-1))/(1+exp(r.n3))

  # variance of expression levels for cluster 3 for normal subjects
  delta.n3<-paraRP2[13]; 
  s.n3<-s.c3+exp(delta.n3)
  sigma2.n3<-exp(s.n3)

  ##############################
  para<- c(pi.1, pi.2, pi.3,
              sigma2.c1, sigma2.n1, mu.c1, rho.c1, mu.n1, rho.n1,
              sigma2.2, mu.c2, rho.c2, mu.n2, rho.n2,
              sigma2.c3, sigma2.n3, mu.c3, rho.c3, mu.n3, rho.n3)

  names(para)<-paraNames

  return(para)
}


