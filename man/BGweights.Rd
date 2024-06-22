\encoding{utf-8}
\name{BGWeights}
\alias{BGWeights}
\title{Bates-Granger minimal variance model weights}
\usage{
BGWeights(object, ..., data, force.update = FALSE)
}
\arguments{
\item{object, \dots}{two or more fitted \code{\link{glm}} objects, or a
\code{list} of such, or an \code{\link[=model.avg]{"averaging"}} object.}

\item{data}{a data frame containing the variables in the model.}

\item{force.update}{if \code{TRUE}, the much less efficient method of
updating \code{glm} function will be  used rather than directly \emph{via}
\code{\link{glm.fit}}. This only applies to \code{glm}s, in
case of other model types \code{update} is always used.}
}
\value{
A numeric vector of model weights.
}
\description{
Compute empirical weights based on out of sample forecast variances,
following Bates and Granger (1969).
}
\details{
Bates-Granger model weights are calculated using prediction covariance. To
get the estimate of prediction covariance, the models are fitted to
randomly selected half of \code{data} and prediction is done on the
remaining half.
These predictions are then used to compute the variance-covariance between
models, \eqn{\Sigma}. Model weights are then calculated as
\myequation{w_{BG} = (1' \Sigma^{-1} 1)^{-1} 1 \Sigma^{-1}
}{w_BG = (1' \Sigma{^-1} 1){^-1} 1 \Sigma{^-1}
}{w_BG = (1' Sigma^-1 1)^-1 1 \ Sigma^-1},
where \eqn{1} a vector of 1-s.

Bates-Granger model weights may be outside of the \eqn{[0,1]} range, which
may cause the averaged variances to be negative. Apparently this method
works best when data is large.
}
\note{
For matrix inversion, \code{\link[MASS:ginv]{MASS::ginv()}} is more stable near singularities
than \code{\link{solve}}. It will be used as a fallback if \code{solve} fails and
\pkg{MASS} is available.
}
\examples{
fm <- glm(Prop ~ mortality + dose, family = binomial, Beetle, na.action = na.fail)
models <- lapply(dredge(fm, evaluate = FALSE), eval)
ma <- model.avg(models)

# this produces warnings because of negative variances:
set.seed(78)
Weights(ma) <- BGWeights(ma, data = Beetle)
coefTable(ma, full = TRUE)

# SE for prediction is not reliable if some or none of coefficient's SE
# are available
predict(ma, data = test.data, se.fit = TRUE)
coefTable(ma, full = TRUE)

}
\references{
Bates, J. M. and Granger, C. W. J. 1969 The combination of forecasts.
\emph{Journal of the Operational Research Society} \strong{20}, 451-468.

Dormann, C. et al. (2018) Model averaging in ecology: a review of Bayesian,
information-theoretic, and tactical approaches for predictive inference.
\emph{Ecological Monographs} \strong{88}, 485–504.
}
\seealso{
\code{\link{Weights}}, \code{\link{model.avg}}

Other model weights: 
\code{\link{bootWeights}()},
\code{\link{cos2Weights}()},
\code{\link{jackknifeWeights}()},
\code{\link{stackingWeights}()}
}
\author{
Carsten Dormann, Kamil Barto\enc{ń}{n}
}
\concept{model weights}
\keyword{models}
