\name{GFL}
\alias{GFL}
\alias{Generalized Fused Lasso}

\title{
Solve Generalized Fused Lasso Model
}

\description{
This function use MM algorithm to fit a piece-wise constant curve to each of the signal sequences.
}
\usage{
GFL(Y, Delta, sigma, rho1 = 1, rho2 = 2, rho3 = 0, 
    obj_c = NULL, max_iter = 1000, verbose = FALSE)
}

\arguments{
  \item{Y}{
    A matrix of original signal, where each column corresponds to a sequence and each row correspond to a marker.
  }
  \item{Delta}{
    A matrix of the same dimension as \code{Y}. Each entry is associated with the entry of \code{Y} at the same location, 
    indicating whether the value of the corresponding entry is regarded as missing or not. 1 = available, 0 = missing. 
  }
  \item{sigma}{
    A vector of standard deviations for each column of \code{Y}. It should have the same length as the number of columns of \code{Y}.
  }
  \item{rho1,rho2,rho3}{
    Factors to be set in the tuning parameters of lambda1, lambda2, and lambda3. See details.
  }
  \item{obj_c}{
    Stopping criterion based on the size of improvement of objective function.
  }
  \item{max_iter}{
    Maximum iteration of MM algorithm to be used to solve the GFL model.
  }
  \item{verbose}{
    Logical. It indicates whether display the intermediate diagnosis imformation. Defautl is \code{FALSE} (highly recommended).
  }
}

\details{
In order to fit a piece-wise constant curve to each of the signal sequences, we try to minimize the following objective function
  \deqn{\frac{1}{2}\sum_{i=1}^M\sum_{j=1}^N (y_{ij}\delta_{ij}-\beta_{ij}\delta_{ij})^2 + \sum_{i=1}^M \lambda_{1,i} \sum_{j=1}^N |\beta_{ij}|
        +\sum_{i=1}^M\lambda_{2,i}\sum_{j=2}^N |\beta_{ij} - \beta_{i,j-1}| 
        +\sum_{j=2}^N \left[ \sum_{i=1}^M \lambda_{3,i}^2 (\beta_{ij}-\beta_{i,j-1})^2 \right]^{\frac{1}{2}}}{loss function
        + lambda1 * lasso penalty + lambda2 * fused lasso penalty + lambda3 * group fused lasso penalty
    }
The optimal solution is approached via an iteration based algorithm called Majorization-Minimization (MM) algorithm developed by Kenneth Lange (2004). The choices of tuning parameters of the model are suggested as follows:
  \deqn{\lambda_{1,i} = c_1 \sigma_i}
  \deqn{\lambda_{2,i} = \rho(p) c_2 \sigma_i \sqrt{\log N}}
  \deqn{\lambda_{3,i} = [1-\rho(p)] c_3 \sigma_i \sqrt{pM} \sqrt{\log N}}
where \eqn{\sigma_i} is signal noise level of each sequence, \eqn{M} is the number of sequences, \eqn{N} is the number of markers and \eqn{c_1}, \eqn{c_2}, \eqn{c_3}, \eqn{\rho}, and \eqn{p} are properly chosen contants, which are absorbed in \eqn{\rho_1}, \eqn{\rho_2}, and \eqn{\rho_3} respectively. More details are referred to Zhang et al. (2012).
}

\value{
  All outputs are collected in a list: 
  \item{obj}{A vector of values of objective function at each MM iteration.}
  \item{Beta}{A matrix of the same dimension as \code{Y}, recording the fitted piece-wise contant curves for each sequence. One column
  	correpond to one sequence, while one row reprents one marker.}
}

\references{
  \enumerate{
  	\item Kenneth Lange. (2004) \emph{Optimization}. Springer, New York.
    \item Zhongyang Zhang, Kenneth Lange, Roel Ophoff, and Chiara Sabatti. (2010) Reconstructing DNA copy number by penalized estimation and imputation. \emph{The Annals of Applied Statistics}, 4(4): 1749-1773.
    \item Zhongyang Zhang, Kenneth Lange, and Chiara Sabatti. (2012) Reconstructing DNA copy number by segmentation of multiple sequences. \emph{Submitted}.
  }
}

\author{
    Zhongyang (Thomas) Zhang, \email{zhangzy@ucla.edu}
    
    %% Kenneth Lang, \email{klange@ucla.edu}
    
    %% Chiara Sabatti, \email{sabatti@stanford.edu}
}

\note{
\code{Y} and \code{Delta} must be of the class matrix. If only one signal sequence is to be analyzed, they should be also coerced to matrix with only one column.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{\link{FL}} for segmentation of only one sequence of signals.
}

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

\keyword{Fused Lasso}
\keyword{MM algorithm}
