\name{BAF.transform}
\alias{BAF.transform}
\alias{mBAF}

\title{
Transform BAF into mBAF
}

\description{
This function is dedicated to transform BAF value into mirrored BAF (mBAF) value. Non-informative SNPs for 
CNV inference have been removed, while missing values for those removed SNPs are initialized with 
the average of nearest SNPs.
}

\usage{
BAF.transform(x, gt = NULL, mBAF.thd = 0.97, win.thd = 0.8, 
              w = 1, k = 2, median.adjust = FALSE)
}

\arguments{
  \item{x}{
    A vector of BAF values to be transformed.
  }
  \item{gt}{
    In tumor data set, if the tumor sample under investigation has matched normal tissue sample, 
    \code{gt} indicates the vector of the genotypes of SNPs in matched normal sample. If no such 
    information can be supplied, it is set \code{NULL} as default.
  }
  \item{mBAF.thd}{
    A criteria to remove non-informative SNPs if no information from matched normal tissue is supplied. 
    See reference for more details.
  }
  \item{win.thd}{
    A further criteria to remove possible non-informative SNPs which might pass the \code{mBAF.thd} criteria. 
    See reference for more details.
  }
  \item{w}{
    The window size used in computation of a quantity to be compared with \code{win.thd}. The default is \eqn{1}. 
    See reference for more details.
  }
  \item{k}{
    The number of nearest SNPs used to computed the initialized values of removed non-informative SNPs.
  }
  \item{median.adjust}{
    Logical. If it is \code{TRUE}, the median of BAF value in between \eqn{0.25} and \eqn{0.75} will be adjusted to 
    \eqn{0.5} first before any transformation applied.
  }
}

\details{
  More details about the transformation are referred to Staaf J., et al. (2008). The missing values for removed non-informative 
  SNPs are initialized with the average of k-nearest SNPs plus a normal random noise in order to eliminate the dependence of adjacent
  SNPs.
}

\value{
  All returned information is collected into a list
  \item{mBAF}{A vector of mirrored BAF values. Missing values of removed non-informative SNPs are initialized for downstream analysis.}
  \item{idx}{A vector of indices of those informative SNPs with values remaining after transformation.}
  \item{idx.na}{A vector of indices of those non-informative SNPs with orignal values removed.}
}

\references{
Staaf J.,  et al. (2008) Segmentation-based detection of allelic imbalance and loss-of-heterozygosity
in cancer cells using whole genome SNP arrays. \emph{Genome Biology}, 9: R136+.
}

\author{
    Zhongyang (Thomas) Zhang, \email{zhangzy@ucla.edu}
    
%%    Kenneth Lang, \email{klange@ucla.edu}
    
%%    Chiara Sabatti, \email{sabatti@stanford.edu}
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
## simulate a sequence of BAF values for 100 SNPs
xf <- sample(x=c(0,0.5,1),size=100,replace=TRUE,prob=c(0.25,0.5,0.25)) + rnorm(100,0,0.02)
xf[xf<0] <- 0
xf[xf>1] <- 1
## insert the signal pattern of a duplcation in the middle of x1
xm <- sample(x=c(0,1),size=20,replace=TRUE,prob=c(0.5,0.5)) + rnorm(20,0,0.02)
xm[xm<0] <- 0
xm[xm>1] <- 1
xf[41:60] <- 2/3*xf[41:60] + 1/3*xm
BAF <- xf
plot(BAF,xlab="SNP",ylab="BAF")

## tranform BAF to mBAF
res <- BAF.transform(x=BAF, gt = NULL, mBAF.thd = 0.97, win.thd = 0.8, 
              w = 1, k = 2, median.adjust = FALSE)
plot(res$mBAF,type="n",xlab="SNP",ylab="mBAF")
points(res$idx,res$mBAF[res$idx])
points(res$idx.na,res$mBAF[res$idx.na],col="gray")
}

\keyword{BAF}
\keyword{mBAF}
