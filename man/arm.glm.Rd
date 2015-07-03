\name{arm.glm}
\alias{arm.glm}
\encoding{utf-8}
\title{Adaptive Regression by Mixing}
\description{
Combine all-subsets GLMs using the ARM algorithm.
}

\usage{
arm.glm(object, R = 250, weight.by = c("aic", "loglik"), trace = FALSE)
}

\arguments{
  \item{object}{a fitted \dQuote{global} \code{glm} object. }
  \item{R}{number of permutations. }
  \item{weight.by}{indicates whether model weights should be calculated with AIC
    or log-likelihood. }
  \item{trace}{if \code{TRUE}, information is printed during the running of
    \code{arm.glm}. }
}

\details{
For each of all-subsets of the \dQuote{global} model, parameters are estimated 
using randomly sampled half of the data. Log-likelihood given the remaining half
of the data is used to calculate AIC weights. This is repeated \code{R} 
times and mean of the weights is used to average all-subsets parameters 
estimated using complete data.
}

\note{
Number of parameters is limited to \code{floor(nobs(object) / 2) - 1}.
All-subsets respect marginality constraints.
}

\value{
An object of class \code{"averaging"} contaning only \dQuote{full} averaged 
coefficients. See \code{\link{model.avg}} for object description.
}

\references{
Yang Y. (2001) Adaptive Regression by Mixing. 
\emph{Journal of the American Statistical Association} 96: 574–588.

Yang Y. (2003) Regression with multiple candidate models: selecting or mixing? 
\emph{Statistica Sinica} 13: 783–810.
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{model.avg}}, \code{\link{par.avg}}

Other implementation: \code{arms} in (archived) package \bold{MMIX}.
}

\examples{
fm <- glm(y ~ X1 + X2 + X3 + X4, data = Cement)

summary(arm.glm(fm, R = 25))
}

\keyword{models}
