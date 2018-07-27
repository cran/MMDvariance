\name{plotHistDensity.v}
\alias{plotHistDensity.v}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot of histogram and density estimate of the pooled gene expression levels. }
\description{
Plot of histogram of pooled gene expression levels, composited with density estimate based on the mixture of marginal distributions. The density estimate is based on the assumption that the marginal correlations between subjects are zero.
}
\usage{
plotHistDensity.v(obj.gsMMD,
                plotFlag="case",
                plotComponent=FALSE,
                myxlab="expression level",
                myylab="density",
                mytitle="Histogram (case)",
                x.legend=NULL,
                y.legend=NULL,
                numPoints=500,
                mycol=1:4, 
                mylty=1:4, 
                mylwd=rep(3,4), 
                cex.main=2, 
                cex.lab=1.5, 
                cex.axis=1.5, 
                cex=2,
                bty="n")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.gsMMD}{an object returned by \code{gsMMD.v},
\code{gsMMD.default.v}, \code{gsMMD2.v},
or \code{gsMMD2.default.v}
}
  \item{plotFlag}{logical. Indicate the plot will based on which type of subjects.}
  \item{plotComponent}{logical. Indicate if components of the mixture of marginal distribution will be plotted.}
  \item{myxlab}{label for x-axis}
  \item{myylab}{label for y-axis}
  \item{mytitle}{title of the plot}
  \item{x.legend}{the x-corrdiates of the legend}
  \item{y.legend}{the y-corrdiates of the legend}
  \item{numPoints}{logical. Indicate how many genes will be plots.}
  \item{mycol}{color for the density estimates (overall and components)}
  \item{mylty}{line styles for the density estimates (overall and components)}
  \item{mylwd}{line width for the density estimates (overall and components)}
  \item{cex.main}{font for main title}
  \item{cex.lab}{font for x- and y-axis labels}
  \item{cex.axis}{font for x- and y-axis}
  \item{cex}{font for texts}
  \item{bty}{the type of box to be drawn around the legend.  The allowed
          values are '"o"' and '"n"' (the default).}
}
\details{
For a given type of subjects, we pool their expression levels together if
the marginal correlations among subjects are zero. We then draw a histogram
of the pooled expression levels. Next, we composite density estimates of gene
expression levels for the overal distribution and the 3 component distributions.
}
\value{
A list containing coordinates of the density estimates:
  \item{x }{sorted pooled gene expression levels for cases or controls.}
  \item{x2 }{a subset of \code{x} specified by the sequence:
\code{seq(from=1,to=len.x, by=delta)}, where \code{len.x} is the
length of the vector \code{x}, and \code{delta=floor(len.x/numpoints).}}
  \item{y }{density estimate corresponding to \code{x2}}
  \item{y1 }{weighted density estimate for gene cluster 1}
  \item{y2 }{weighted density estimate for gene cluster 2}
  \item{y3 }{weighted density estimate for gene cluster 3}
}
\references{
  Li X, Fu Y, Wang X, DeMeo DL, Tantisira K, Weiss ST, Qiu W.
  Detecting Differentially Variable MicroRNAs via Model-Based Clustering.
  \emph{International Journal of Genomics}. Article ID 6591634, Volumne 2018 (2018).
}
\author{ 
Xuan Li \email{lixuan0759@gmail.com},
Yuejiao Fu \email{yuejiao@mathstat.yorku.ca},
Xiaogang Wang \email{stevenw@mathstat.yorku.ca},
Dawn L. DeMeo \email{redld@channing.harvard.edu},
Kelan Tantisira \email{rekgt@channing.harvard.edu},
Scott T. Weiss \email{restw@channing.harvard.edu},
Weiliang Qiu \email{weiliang.qiu@gmail.com}
}
\note{ The density estimate is obtained based on the
assumption that the marginal correlation among
subjects is zero. If the estimated marginal correlation obtained by \code{gsMMD.v} is far from zero, then do not use this plot function. 
}
\examples{
  
      t1 = proc.time()
      library(ALL)
      data(ALL)
      eSet1 <- ALL[1:50, ALL$BT == "B3" | ALL$BT == "T2"]
      
      mem.str <- as.character(eSet1$BT)
      nSubjects <- length(mem.str)
      memSubjects <- rep(0,nSubjects)
      # B3 coded as 0, T2 coded as 1
      memSubjects[mem.str == "T2"] <- 1
      
      obj.gsMMD.v <- gsMMD.v(eSet1, memSubjects, transformFlag = FALSE, 
        transformMethod = "boxcox", scaleFlag = FALSE, 
        eps = 1.0e-1, ITMAX = 5, quiet = TRUE)
      print(round(obj.gsMMD.v$para, 3))
     
  
    plotHistDensity.v(obj.gsMMD.v, plotFlag = "case", 
        mytitle = "Histogram (case)", 
        plotComponent = TRUE, 
        x.legend = c(0.8, 3), 
        y.legend = c(0.3, 0.4), 
        numPoints = 50)
    t2=proc.time()-t1
    print(t2)
  

}
\keyword{classif }% at least one, from doc/KEYWORDS
