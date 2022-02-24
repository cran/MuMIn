\encoding{utf-8}
\name{stackingWeights}
\alias{stackingWeights}
\title{Stacking model weights}
\usage{
stackingWeights(object, ..., data, R, p = 0.5)
}
\arguments{
\item{object, \dots}{two or more fitted \code{\link{glm}} objects, or a
\code{list} of such, or an \code{\link[=model.avg]{"averaging"}} object.}

\item{data}{a data frame containing the variables in the model, used for
fitting and prediction.}

\item{R}{the number of replicates.}

\item{p}{the proportion of the \code{data} to be used as training set.
Defaults to 0.5.}
}
\value{
A matrix with two rows, containing model weights
calculated using \code{mean} and \code{median}.
}
\description{
Computes model weights based on a cross-validation-like procedure.
}
\details{
Each model in a set is fitted to the training data: a subset of \code{p * N}
observations in \code{data}. From these models a prediction is produced on
the remaining part of \code{data} (the test
or hold-out data). These hold-out predictions are fitted to the hold-out
observations, by optimising the weights by which the models are combined. This
process is repeated \code{R} times, yielding a distribution of weights for each
model (which Smyth & Wolpert (1998) referred to as an \sQuote{empirical Bayesian
estimate of posterior model probability}). A mean or median of model weights for
each model is taken and re-scaled to sum to one.
}
\note{
This approach requires a sample size of at least \eqn{2\times}{2x} the number
of models.
}
\examples{
#simulated Cement dataset to increase sample size for the training data 
fm0 <- glm(y ~ X1 + X2 + X3 + X4, data = Cement, na.action = na.fail)
dat <- as.data.frame(apply(Cement[, -1], 2, sample, 50, replace = TRUE))
dat$y <- rnorm(nrow(dat), predict(fm0), sigma(fm0))

# global model fitted to training data:
fm <- glm(y ~ X1 + X2 + X3 + X4, data = dat, na.action = na.fail)

# generate a list of *some* subsets of the global model
models <- lapply(dredge(fm, evaluate = FALSE, fixed = "X1", m.lim = c(1, 3)), eval)

wts <- stackingWeights(models, data = dat, R = 10)

ma <- model.avg(models)
Weights(ma) <- wts["mean", ]

predict(ma)

}
\references{
Wolpert, D. H. (1992) Stacked generalization. \emph{Neural Networks}, 5: 241-259.

Smyth, P. & Wolpert, D. (1998) \emph{An Evaluation of Linearly Combining
Density Estimators via Stacking. Technical Report No. 98-25.} Information
and Computer Science Department, University of California, Irvine, CA.

Dormann, C. et al. (2018) Model averaging in ecology: a review of Bayesian,
information-theoretic, and tactical approaches for predictive inference.
\emph{Ecological Monographs}, 88, 485–504.
}
\seealso{
\code{\link{Weights}}, \code{\link{model.avg}}

Other model weights: 
\code{\link{BGWeights}()},
\code{\link{bootWeights}()},
\code{\link{cos2Weights}()},
\code{\link{jackknifeWeights}()}
}
\author{
Carsten Dormann, Kamil Barto\enc{ń}{n}
}
\concept{model weights}
\keyword{models}
