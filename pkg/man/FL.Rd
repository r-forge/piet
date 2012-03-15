\name{FL}
\alias{FL}
\alias{Fused Lasso}

\title{
Solve Fused Lasso Model
}
\description{
This function uses MM algorithm to fit a piece-wise constant curve to a sequence of signals ordered linearly.
}

\usage{
FL(y, sigma, rho1 = 1, rho2 = 2, obj_c = 1e-04, max_iter = 1000)
}


\arguments{
  \item{y}{
    A vector of original signal.
  }
  \item{sigma}{
    A numeric number indicating the standard deviation of \code{y}
  }
  \item{rho1,rho2}{
    Factors to be set in the tuning parameters of lambda1 and lambda2. See details.
  }
  \item{obj_c}{
    Stopping criterion based on the size of improvement of objective function.
  }
  \item{max_iter}{
    Maximum iteration of MM algorithm to be used to solve the GFL model.
  }
}

\details{
In order to fit a piece-wise constant curve to signal intensities ordered linearly, we try to minimize the following objective function
  \deqn{\frac{1}{2}\sum_{i=1}^N (y_i-\beta_i)^2 + \lambda_1 \sum_{i=1}^N |\beta_i|
        +\lambda_2\sum_{i=2}^N |\beta_i - \beta_{i-1}| }{
        loss function + lambda1 * lasso penalty + lambda2 * fused lasso penalty
    }
The optimal solution is approached via an iteration based algorithm called Majorization-Minimization (MM) algorithm developed by Kenneth Lange (2004).The choices of tuning parameters of the model are suggested as follows:
  \deqn{\lambda_1 = \rho_1 \sigma}
  \deqn{\lambda_2 = \rho_2 \sigma \sqrt{\log N}}
where \eqn{\sigma} is an estimate of standard deviation of signals, \eqn{N} is the number of markers and \eqn{\rho_1} and \eqn{\rho_2} are properly chosen contant factors. More details are referred to Zhang et al. (2010).
}

\value{
  All outputs are collected in a list: 
  \item{obj}{A vector of values of objective function at each MM iteration.}
  \item{beta}{A vector of the same dimension as \code{y}, recording the fitted piece-wise contant curve.}
}

\references{
  \enumerate{
  	\item Kenneth Lange. (2004) \emph{Optimization}. Springer, New York.
    \item Zhongyang Zhang, Kenneth Lange, Roel Ophoff, and Chiara Sabatti. (2010) Reconstructing DNA copy number by penalized estimation and imputation. \emph{The Annals of Applied Statistics}, 4(4): 1749-1773.
  }
}

\author{
    Zhongyang (Thomas) Zhang, \email{zhangzy@ucla.edu}

%%    Kenneth Lang, \email{klange@ucla.edu}
    
%%    Chiara Sabatti, \email{sabatti@stanford.edu}
}

\note{
If the user just wants to segment one sequence of signals, then \code{GL} is a little more efficient than \code{\link{GFL}}, which is designed for segmentation of multiple sequences. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{GFL}} for joint segmentation of multiple sequences of signals.
}

\examples{
## Segment 1 sequence of signals with 100 markers
## Duplications are superimposed in the middle
y <- rnorm(100,0,0.15)
y[41:60] <- rnorm(20,0.3,0.2)
tmp <- y[-1] - y[-length(y)]
sigma <- sd(tmp)/sqrt(2)
res <- FL(y=y, sigma=sigma, rho1 = 1, rho2 = 2, obj_c = 1e-04, max_iter = 1000)
plot(res$beta,type="s")
}

\keyword{Fused Lasso}
\keyword{MM algorithm}
