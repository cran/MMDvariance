\name{gsMMD.v}
\alias{gsMMD.v}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gene selection based on variances by using a mixture of marginal distributions}
\description{
Gene selection based on variances by using the marginal distributions of gene profiles that characterized by a mixture of three-component multivariate distributions. The goal is to detect gene probes having different variances between
cases and controls.
Input is an object derived from the class \code{ExpressionSet}. The function will obtain initial gene cluster membership by its own.
}
\usage{
gsMMD.v(obj.eSet, 
      memSubjects, 
      maxFlag = TRUE, 
      thrshPostProb = 0.5, 
      geneNames = NULL, 
      alpha = 0.05, 
      iniGeneMethod = "myLeveneTest", 
      transformFlag = FALSE, 
      transformMethod = "boxcox", 
      scaleFlag = TRUE, 
      criterion = c("cor", "skewness", "kurtosis"), 
      minL = -10, 
      maxL = 10, 
      stepL = 0.1, 
      eps = 0.001, 
      ITMAX = 100, 
      plotFlag = FALSE,
      quiet=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{obj.eSet}{an object derived from the class \code{ExpressionSet} which contains the matrix of gene expression levels.  The rows of the matrix are genes. The columns of the matrix are subjects.}
  \item{memSubjects}{a vector of membership of subjects. \code{memSubjects[i]=1} means the \eqn{i}-th subject belongs to diseased group, \eqn{0} otherwise. 
  }
  \item{maxFlag}{logical. Indicate how to assign gene class membership. \code{maxFlag}=TRUE means that a gene will be assigned
to a class in which the posterior probability of the gene belongs to this class is maximum. \code{maxFlag}=FALSE means that a gene will be assigned to class 1 if the posterior probability of the gene belongs to class 1 is
greater than \code{thrshPostProb}. Similarly, a gene will be assigned to class 1 if the posterior probability of the gene belongs to class 1 is
greater than \code{thrshPostProb}. If  the posterior probability is less than \code{thrshPostProb}, the gene
will be assigned to class 2 (non-differentially variable gene group).}
  \item{thrshPostProb}{threshold for posterior probabilities. For example, if the posterior probability that a gene belongs to cluster 1 given its gene expression levels is larger than \code{thrshPostProb}, then this gene will be assigned to cluster 1.}
  \item{geneNames}{an optional character vector of gene names}
  \item{alpha}{significant level which is equal to \code{1-conf.level}, 
\code{conf.level} is the argument for the function \code{t.test}.
  }
  \item{iniGeneMethod}{method to get initial 3-cluster partition of genes: 
(1) genes having higher variance in cases than in controls;
(2) genes having equal variance between cases and controls;
(3) genes having lower variance in cases than in controls.

Available methods are: 
\dQuote{myAWvar}, 
\dQuote{myBFTest}, 
\dQuote{myFTest}, 
\dQuote{myLeveneTest}, 
\dQuote{myLevene.TM}, 
\dQuote{myiAWvar.BF},
\dQuote{myiAWvar.Levene}, 
\dQuote{myiAWvar.TM},
\dQuote{myLeveneTest}, 
\dQuote{myLeveneTest.TM}. 
  }
  \item{transformFlag}{logical. Indicate if data transformation is needed}
  \item{transformMethod}{method for transforming data. Available methods include "boxcox", "log2", "log10", "log", "none".}
  \item{scaleFlag}{logical. Indicate if gene profiles are to be scaled to have mean zero and variance one. If \code{transformFlag=TRUE} and \code{scaleFlag=TRUE}, then scaling is performed after transformation. To avoid linear dependence of tissue samples after scaling
    gene profiles, we delete one tissue sample after scaling (c.f. details).}
  \item{criterion}{if \code{transformFlag=TRUE}, \code{criterion} indicates what criterion to determine if data looks like normal.  \dQuote{cor} means using Pearson's correlation. The idea is that the observed quantiles after transformation should be close to theoretical normal quantiles. So we can use Pearson's correlation to check if the scatter plot of theoretical normal quantiles versus observed quantiles is a straightline.  \dQuote{skewness} means using skewness measure to check if the distribution of the transformed data are close to normal distribution; \dQuote{kurtosis} means using kurtosis measure to check normality.}
  \item{minL}{lower limit for the \code{lambda} parameter used in Box-Cox transformation}
  \item{maxL}{upper limit for the \code{lambda} parameter used in Box-Cox transformation}
  \item{stepL}{step increase when searching the optimal \code{lambda} parameter used in Box-Cox transformation}
  \item{eps}{a small positive value. If the absolute value of a value is smaller than \code{eps}, this value is regarded as zero.  }
  \item{ITMAX}{maximum iteration allowed for iterations in the EM algorithm}
  \item{plotFlag}{logical. Indicate if the Box-Cox normality plot should be output.}
  \item{quiet}{logical. Indicate if intermediate results should be printed out.}
}
\details{
We assume that the distribution of gene expression profiles is 
a mixture of 3-component multivariate normal distributions 
\eqn{\sum_{k=1}^{3} \pi_k f_k(x|\theta)}. Each component distribution \eqn{f_k} 
corresponds to a gene cluster. The 3 components correspond to 3 gene clusters:
(1) genes having higher variance in cases than in controls;
(2) genes having equal variance between cases and controls;
(3) genes having lower variance in cases than in controls.
The model parameter vector is \eqn{\theta=(\pi_1}, \eqn{\pi_2}, 
\eqn{\pi_3}, 
\eqn{\sigma^2_{c1}}, \eqn{\sigma^2_{n1}}, \eqn{\mu_{c1}}, \eqn{\rho_{c1}}, \eqn{\mu_{n1}}, 
\eqn{\rho_{n1}},
\eqn{\sigma^2_2}, \eqn{\mu_{c2}}, \eqn{\rho_{c2}}, \eqn{\mu_{n2}}, \eqn{\rho_{n2}}, 
\eqn{\sigma^2_{c3}}, \eqn{\sigma^2_{n3}}, 
\eqn{\mu_{c3}}, 
\eqn{\rho_{c3}}, \eqn{\mu_{n3}}, \eqn{\rho_{n3}}.
where \eqn{\pi_1}, \eqn{\pi_2}, and \eqn{\pi_3} are the mixing proportions; 
\eqn{\mu_{c1}}, \eqn{\sigma^2_{c1}}, and \eqn{\rho_{c1}} are 
the marginal mean, variance, and correlation of gene expression levels 
of cluster 1 (over-variable genes) for diseased subjects; 
\eqn{\mu_{n1}}, \eqn{\sigma^2_{n1}}, and \eqn{\rho_{n1}} are 
the marginal mean, variance, and correlation of gene expression levels 
of cluster 1 (over-variable genes) for non-diseased subjects; 
\eqn{\sigma^2_2}, \eqn{\mu_{c2}}, \eqn{\rho_{c2}}, \eqn{\mu_{n2}}, and \eqn{\rho_{n2}} are the marginal mean, 
variance, and correlation of gene 
expression levels of cluster 2 (equal-variable genes); 
\eqn{\mu_{c3}}, \eqn{\sigma^2_{c3}}, and \eqn{\rho_{c3}} are 
the marginal mean, variance, and correlation of gene expression levels 
of cluster 3 (under-variable genes) for diseased subjects; 
\eqn{\mu_{n3}}, \eqn{\sigma^2_{n3}}, and \eqn{\rho_{n3}} are 
the marginal mean, variance, and correlation of gene expression levels 
of cluster 3 (under-variable) for non-diseased subjects. 

Note that genes in cluster 2 are non-differentially variable across
abnormal and normal tissue samples. Hence there are only 5 parameters for
cluster 2.

To make sure the identifiability, we set the following contraints: 
\eqn{\sigma_{c1}>\sigma_{n1}} and \eqn{\sigma_{c3}<\sigma_{n3}}.

To make sure the marginal covariance matrices are poisitive definite,
we set the following contraints: 
\eqn{-1/(n_c-1)<\rho_{c1}<1},
\eqn{-1/(n_n-1)<\rho_{n1}<1},
\eqn{-1/(n-1)<\rho_{2}<1},
\eqn{-1/(n_c-1)<\rho_{c3}<1},
\eqn{-1/(n_n-1)<\rho_{n3}<1}.

We also has the following constraints for the mixing proportion:
\eqn{\pi_3=1-\pi_1-\pi_2}, \eqn{\pi_k>0}, \eqn{k=1,2,3}.

We apply the EM algorithm to estimate the model parameters. 
We regard the cluster membership of genes as missing values.

To facilitate the estimation of the parameters,
we reparametrize the parameter vector as
\eqn{\theta^*=(\pi_1}, \eqn{\pi_2}, 
\eqn{s^2_{c1}}, \eqn{\delta_{n1}}, \eqn{\mu_{c1}}, 
\eqn{r_{c1}}, \eqn{\mu_{n1}}, \eqn{r_{n1}},
\eqn{s^2_2}, \eqn{\mu_{c2}}, \eqn{r_{c2}}, \eqn{\mu_{n2}}, \eqn{r_{n2}}, 
\eqn{s^2_{c3}}, \eqn{\delta_{n3}}, \eqn{\mu_{c3}}, \eqn{r_{c3}}, \eqn{\mu_{n3}}, 
\eqn{r_{n3})},
where
\eqn{\sigma_{n1}=\sigma_{c1}-\exp(\delta_{n1})},
\eqn{\sigma_{n3}=\sigma_{c3}+\exp(\delta_{n3})},
\eqn{\rho_{c1}=(\exp(r_{c1})-1/(n_c-1))/(1+\exp(r_{c1}))},
\eqn{\rho_{n1}=(\exp(r_{n1})-1/(n_n-1))/(1+\exp(r_{n1}))},
\eqn{\rho_{2}=(\exp(r_{2})-1/(n-1))/(1+\exp(r_{2}))},
\eqn{\rho_{c3}=(\exp(r_{c3})-1/(n_c-1))/(1+\exp(r_{c3}))},
\eqn{\rho_{n3}=(\exp(r_{n3})-1/(n_n-1))/(1+\exp(r_{n3}))}.

Given a gene, the expression levels of the gene are assumed independent. However, after scaling, the scaled expression levels of the gene are no longer independent and the rank \eqn{r^*=r-1} of the covariance matrix for the scaled gene profile will be one less than the rank \eqn{r} for the un-scaled gene profile Hence the covariance matrix of the
gene profile will no longer be positive-definite. 
To avoid this problem,
we delete a tissue sample after scaling since its information has been incorrporated by other scaled tissue samples. We arbitrarily select the tissue sample, which has the biggest label number, from the tissue sample group that has larger size than the other tissue sample group. For example, if there are 6 cancer tissue samples and 10 normal tissue samples, we delete the 10-th normal tissue sample after scaling.
}
\value{
A list contains 18 elements.
  \item{dat}{the (transformed) microarray data matrix. If tranformation
performed, then \code{dat} will be different from the input 
microarray data matrix.}
  \item{memSubjects}{the same as the input \code{memSubjects}.}
  \item{memGenes}{a vector of cluster membership of genes. \eqn{1} means over-variable gene; \eqn{2} means non-differentially variable gene; \eqn{3} 
means under-variable gene.}
  \item{memGenes2}{an variant of the vector of cluster membership of genes. 
\eqn{1} means differentially variable gene; \eqn{0} means non-differentially variable gene.}
  \item{para}{parameter estimates (c.f. details).}
  \item{llkh}{value of the loglikelihood function.}
  \item{wiMat}{posterior probability that a gene belongs to a cluster given the expression levels of this gene. Column i is for cluster i.}
  \item{wiArray}{posterior probability matrix for different initial gene selection methods.}
  \item{memIniMat}{a matrix of initial cluster membership of genes.}
  \item{paraIniMat}{a matrix of parameter estimates based on initial gene cluster membership.}
  \item{llkhIniVec}{a vector of values of loglikelihood function.}
  \item{memMat}{a matrix of cluster membership of genes based on the mixture
of marginal models with initial parameter estimates obtained initial gene cluster membership.
  }
  \item{paraMat}{a matrix of parameter estimates based on the mixture
of marginal models with initial parameter estimates obtained initial gene cluster membership.
 }
  \item{llkhVec}{a vector of values of loglikelihood function based on the mixture of marginal models with initial parameter estimates obtained initial gene cluster membership.}
  \item{lambda}{the parameter used to do Box-Cox transformation}
  \item{paraRP}{parameter estimates for reparametrized parameter vector (c.f. details).}
  \item{paraIniMatRP}{a matrix of parameter estimates for reparametrized parameter vector based on initial gene cluster membership.}
  \item{paraMatRP}{a matrix of parameter estimates for reparametrized parameter vector based on the mixture
of marginal models with initial parameter estimates obtained initial gene cluster membership.}
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
%\seealso{
%}
\note{
The speed of the program is slow for large data sets.
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
      t2=proc.time()-t1
      print(t2)
      
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classif }
