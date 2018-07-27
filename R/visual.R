# created on June 28, 2012
#
plotHistDensity.v<-function(obj.gsMMD,
                          plotFlag="case",
                          plotComponent=FALSE,
                          myxlab="expression level",
                          myylab="density",
                          mytitle="Histogram (case)",
                          x.legend=NULL, y.legend=NULL,
                          numPoints=500,
                          mycol=1:4, mylty=1:4, mylwd=rep(3,4),
                          cex.main=2, cex.lab=1.5, cex.axis=1.5, cex=2,
                          bty="n")
{
  plotFlag<-match.arg(plotFlag, choices=c("case", "control"))

  dat<-obj.gsMMD$dat
  para<-as.numeric(obj.gsMMD$para)
  names(para)<-paraNames
  memSubjects<-obj.gsMMD$memSubjects

  nCases<-sum(memSubjects==1, na.rm=TRUE)
  nControls<-sum(memSubjects==0, na.rm=TRUE)
  n<-nCases+nControls

  pi.1<-para[1]
  pi.2<-para[2]
  pi.3<-para[3]

  sigma2.c1<-para[4]
  sigma2.n1<-para[5]
  mu.c1<-para[6]
  rho.c1<-para[7]
  mu.n1<-para[8]
  rho.n1<-para[9]

  sigma2.2<-para[10]
  mu.c2<-para[11]
  rho.c2<-para[12]
  mu.n2<-para[13]
  rho.n2<-para[14]

  sigma2.c3<-para[15]
  sigma2.n3<-para[16]
  mu.c3<-para[17]
  rho.c3<-para[18]
  mu.n3<-para[19]
  rho.n3<-para[20]

  if(plotFlag=="case")
  { x<-as.vector(dat[,memSubjects==1, drop=FALSE])
  }
  else {
    x<-as.vector(dat[,memSubjects==0, drop=FALSE])
  }
  x<-sort(x)
  len.x<-length(x)
  delta<-floor(len.x/numPoints)
  x2<-x[seq(from=1,to=len.x, by=delta)]

  if(plotFlag=="case")
  { y1<-pi.1*dnorm(x2, mean=mu.c1, sd=sqrt(sigma2.c1))
    y2<-pi.2*dnorm(x2, mean=mu.c2, sd=sqrt(sigma2.2))
    y3<-pi.3*dnorm(x2, mean=mu.c3, sd=sqrt(sigma2.c3))
    y<-y1+y2+y3

  } else {
    y1<-pi.1*dnorm(x2, mean=mu.n1, sd=sqrt(sigma2.n1))
    y2<-pi.2*dnorm(x2, mean=mu.n2, sd=sqrt(sigma2.2))
    y3<-pi.3*dnorm(x2, mean=mu.n3, sd=sqrt(sigma2.n3))
    y<-y1+y2+y3
  }

  tmp<-hist(x, plot=FALSE)
  myylim<-range(c(y, tmp$density), na.rm=TRUE)
  hist(x, freq=FALSE, main=mytitle,xlab=myxlab, ylab=myylab, ylim=myylim,
    cex.main=cex.main, cex.lab=cex.lab)
  lines(x2, y, col=mycol[1], lty=mylty[1], lwd=mylwd[1])

  if(plotComponent)
  {
    lines(x2, y1, col=mycol[2], lty=mylty[2], lwd=mylwd[2])
    lines(x2, y2, col=mycol[3], lty=mylty[3], lwd=mylwd[2])
    lines(x2, y3, col=mycol[4], lty=mylty[4], lwd=mylwd[2])
  }
  if(is.null(x.legend))
  {
    x.max<-max(x, na.rm=TRUE)
    d<-max(x, na.rm=TRUE)-min(x, na.rm=TRUE)
    tmp1<-x.max-d/3
    x.legend<-c(tmp1, x.max)
  }
  if(is.null(y.legend))
  {
    y.max<-max(y, na.rm=TRUE)
    d<-max(y, na.rm=TRUE)-min(y, na.rm=TRUE)
    tmp1<-y.max-d/4
    y.legend<-c(tmp1, y.max)

  }

  if(plotComponent)
  { legend(x=x.legend, y=y.legend, legend=c("overall","component1", "component2"
, "component3"),
    lty=mylty, col=mycol, lwd=mylwd, cex=cex, bty=bty)
  }
  invisible(list(x=x,x2=x2, y=y, y1=y1, y2=y2, y3=y3))
}




