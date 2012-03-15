\name{Piet}
\alias{Piet}
\docType{package}
\title{DNA CNV Analysis Tools Based on Generlized Fused Lasso Model}
\description{
This package provides some segemntation tools for analysis of individual or multiple sequences of CNV data. 
}
\details{
\tabular{ll}{
Package: \tab Piet\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2012-02-05\cr
License: \tab GPL (>=2.0)\cr
}

The most importatn function is \code{GFL} that implements the segmentation algorithm for multiple sequences of signals described in Zhang, Lange, and Sabatti (2012). Another function \code{FL} is an efficient implementation of fused lasso for segmentation of indiviudal sequences. Fitted segment means resulted from \code{GFL} or \code{FL} can be used in subsequence CNV analysis. 
}

\author{
Zhongyang (Thomas) Zhang

Maintainer: Zhongyang (Thomas) Zhang, \email{zhangzy@ucla.edu}
}

\references{
  \enumerate{
    \item Zhongyang Zhang, Kenneth Lange, Roel Ophoff, and Chiara Sabatti. (2010) Reconstructing DNA copy number by penalized estimation and imputation. \emph{The Annals of Applied Statistics}, 4(4): 1749-1773.
    \item Zhongyang Zhang, Kenneth Lange, and Chiara Sabatti. (2012) Reconstructing DNA copy number by segmentation of multiple sequences. \emph{Submitted}.
    }
}

\keyword{ CNV, fused lasso, MM algorithm }

%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}

\examples{
## Jointly segment 2 sequences of signals with 100 markers
## Duplications are superimposed on both sequences
Y <- matrix(rnorm(200,0,0.15),100,2)
Y[41:60,] <- rnorm(40,0.3,0.2)
Delta <- matrix(1,100,2)
sigma <- apply(Y,2,FUN="mad")
res <- GFL(Y=Y, Delta, sigma, rho1 = 0.01, rho2 = 0.5*2, rho3 = 0.5*2, 
           obj_c = 1e-4, max_iter = 1000, verbose = FALSE)
plot(res$Beta[,1],type="s")
}
