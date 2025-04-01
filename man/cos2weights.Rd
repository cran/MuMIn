\encoding{utf-8}
\name{cos2Weights}
\alias{cos2Weights}
\title{Cos-squared model weights}
\usage{
cos2Weights(object, ..., data, eps = 1e-06, maxit = 100, predict.args = list())
}
\arguments{
\item{object, \dots}{two or more fitted \lcode{glm} objects, or a
\code{list} of such, or an \lxcode{"averaging"}{=model.avg} object.
Currently only \code{lm} and \code{glm} objects are accepted.}

\item{data}{a test data frame in which to look for variables
for use with \link[=predict]{prediction}. If omitted, the fitted linear
predictors are used.}

\item{eps}{tolerance for determining convergence.}

\item{maxit}{maximum number of iterations.}

\item{predict.args}{optionally, a \code{list} of additional arguments to be
passed to \code{predict}.}
}
\value{
A numeric vector of model weights.
}
\description{
Calculate the cos-squared model weights, following the algorithm
outlined in the appendix to Garthwaite & Mubwandarikwa (2010).
}
\examples{
\dontshow{
if(length(find.package("expm", quiet = TRUE)) == 1) \{
}
fm <- lm(y ~ X1 + X2 + X3 + X4, Cement, na.action = na.fail)
# most efficient way to produce a list of all-subsets models
models <- lapply(dredge(fm, evaluate = FALSE), eval)
ma <- model.avg(models)

test.data <- Cement
Weights(ma) <- cos2Weights(models, data = test.data)
predict(ma, data = test.data)
\dontshow{
\} else message("Need CRAN package 'expm' to run this example")
}
}
\references{
Garthwaite, P. H. and Mubwandarikwa, E. 2010 Selection of weights for
weighted model averaging. \emph{Australian & New Zealand Journal of
Statistics} \strong{52}, 363–382.

Dormann, C. et al. 2018 Model averaging in ecology: a review of Bayesian,
information-theoretic, and tactical approaches for predictive inference.
\emph{Ecological Monographs} \strong{88}, 485–504.
}
\seealso{
\lcode{Weights}, \lcode{model.avg}

Other model weights: 
\code{\link{BGWeights}()},
\code{\link{bootWeights}()},
\code{\link{jackknifeWeights}()},
\code{\link{stackingWeights}()}
}
\author{
Carsten Dormann, adapted by Kamil Barto\enc{ń}{n}
}
\concept{model weights}
\keyword{models}
