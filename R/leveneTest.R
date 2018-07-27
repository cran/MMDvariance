"leveneTest" <-
function(x, memSubjects, alpha = 0.05, ...)
{
  xc<-x[memSubjects==1]
  xn<-x[memSubjects==0]

  vc<-var(xc, na.rm=TRUE)
  vn<-var(xn, na.rm=TRUE)

  res<-levene.test(y=x, group=as.factor(memSubjects),...)

  if(vc > vn)
  {
    res2<-c(res$p.value, res$statistic)
  } else {
    res2<-c(res$p.value, -res$statistic)
  }
  names(res2)<-c("p.value", "statistic")
  return(res2)
}

