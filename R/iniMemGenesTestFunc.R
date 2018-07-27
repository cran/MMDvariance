"iniMemGenesTestFunc.v" <-
function(X, 
         memSubjects, 
         testFun, 
         geneNames=NULL, 
         alpha=0.05, eps=1.0e-6)
{

  nGenes<-nrow(X)
  alpha.old<-alpha
  #mat<-apply(X, 1, testFun, memSubjects=memSubjects, alpha=alpha)
  mat<-apply(X, 1, testFun, memSubjects=memSubjects)
  mat<-t(mat)
  colnames(mat)<-c("pvalue", "stat")

  # significant and have positive test statistics
  posUpp<-which(mat[,1]<alpha & mat[,2]>0)
  
  # significant and have negative test statistics
  posLow<-which(mat[,1]<alpha & mat[,2]<0)

  len1<-length(posUpp)
  len3<-length(posLow)
  len2<-nGenes-len1-len3
  if(len1>1 && len2>1 && len3>1 && len2>=len1 && len2>=len3)
  { 
  } else {
    mat2<-cbind(1:nrow(mat), mat) 
    mat3<-mat2[order(mat2[,2]),]
    pos.p=which(mat3[,3]>0)
    pos.n=which(mat3[,3]<0)
    if(length(pos.p))
    {
      mat3.p<-mat3[mat3[,3]>0,]
    } else {
      mat3.p=mat3[1:20,]
    }
    if(length(pos.n))
    {
      mat3.n<-mat3[mat3[,3]<0,]
    } else {
      nn=nrow(mat3)
      mat3.n=mat3[c((nn-19):nn),]
    }
    #tmpn<-ceiling(nGenes/5)
    #tmpn = min(nrow(mat3.p), nrow(mat3.n))
    tmpn=10
    #tmpn<-min(max(tmpn, 2, na.rm=TRUE), floor(nGenes/3), na.rm=TRUE)
    #tmpn<-min(c(nrow(mat3.p), nrow(mat3.n),  
    #  max(tmpn, 2, na.rm=TRUE), floor(nGenes/3)), na.rm=TRUE)
    posUpp<-mat3.p[1:tmpn,1]
    posLow<-mat3.n[1:tmpn,1]
  }
  
  memGenes<-rep(2, nGenes)
  memGenes[posUpp]<-1
  memGenes[posLow]<-3

  memGenes2<-rep(1, nGenes)
  memGenes2[memGenes==2]<-0

  if(sum(is.null(geneNames), na.rm=TRUE))
  { geneNames<-paste("gene", 1:nGenes, sep="") }

  names(memGenes)<-geneNames
  names(memGenes2)<-geneNames
  p.value<-mat[,1]
  names(p.value)<-geneNames
  statistic<-mat[,2]
  names(statistic)<-geneNames

  # obtain parameter estimate and value of log likelihood function 
  tmp<-getPara.v(X, memSubjects, memGenes, eps=eps)
  para<-tmp$para
  llkh<-tmp$llkh

  res<-list(para=para, llkh=llkh, memGenes=memGenes, 
            memGenes2=memGenes2, p.value=p.value, 
            statistic=statistic)

  return(res) 
}


