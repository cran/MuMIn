\encoding{utf-8}
\name{jackknifeWeights}
\alias{jackknifeWeights}
\title{Jackknifed model weights}
\usage{
jackknifeWeights(object, ..., data, type = c("loglik", "rmse"),
  family = NULL, weights = NULL, optim.method = "BFGS", maxit = 1000,
  optim.args = list(), start = NULL, force.update = FALSE,
  py.matrix = FALSE)
}
\arguments{
\item{object, \dots}{two or more fitted \code{\link{glm}} objects, or a
\code{list} of such, or an \code{\link[=model.avg]{"averaging"}} object.}

\item{data}{a data frame containing the variables in the model. It is
optional if all models are \code{glm}.}

\item{type}{a character string specifying the function to minimize. Either
\code{"rmse"} or \code{"loglik"}.}

\item{family}{used only if \code{type = "loglik"}, a \code{\link{family}} object to be
used for likelihood calculation. Not needed if all models share the same
\code{family} and link function.}

\item{weights}{an optional vector of \sQuote{\link[=weights.glm]{prior
weights}} to be used in the model fitting process. Should be \code{NULL} or
a numeric vector.}

\item{optim.method}{optional, optimisation method, passed to
\code{\link{optim}}.}

\item{maxit}{optional, the maximum number of iterations, passed to
\code{\link{optim}}.}

\item{optim.args}{optional list of other arguments passed to
\code{\link{optim}}.}

\item{start}{starting values for model weights. Numeric of length equal the
number of models.}

\item{force.update}{for \code{glm}, the \code{glm.fit} function is used for
fitting models to the train data, which is much more efficient. Set to
\code{TRUE} to use \code{update} instead.}

\item{py.matrix}{either a boolean value, then if \code{TRUE} a jackknifed prediction
matrix is returned and if \code{FALSE} a vector of jackknifed model weights, or
a \var{N}×\var{M} matrix (\emph{number of cases} × \emph{number of models})
that is interpreted as a jackknifed prediction matrix and it is used for
optimisation (i.e. the jackknife procedure is skipped).}
}
\value{
The function returns a numeric vector of model weights.
}
\description{
Computes model weights optimized for jackknifed model fits.
}
\details{
Model weights are chosen (using \code{\link{optim}}) to minimise
    \acronym{RMSE} or log-likelihood of
    the prediction for data point \var{i}, of a model fitted omitting that
    data point \var{i}. The jackknife procedure is therefore run for all
    provided models and for all data points.
}
\note{
This procedure can give variable results depending on the 
\link[=optim]{optimisation method} and starting values. It is therefore 
advisable to make several replicates using different \code{optim.method}s. 
See \code{\link{optim}} for possible values for this argument.
}
\examples{
fm <- glm(Prop ~ mortality * dose, binomial(), Beetle, na.action = na.fail)

fits <- lapply(dredge(fm, eval = FALSE), eval)

amJk <- amAICc <- model.avg(fits)
set.seed(666)
Weights(amJk) <- jackknifeWeights(fits, data = Beetle)

coef(amJk)
coef(amAICc)

}
\author{
Kamil Barto\enc{ń}{n}. Carsten Dormann
}
\references{
Hansen, B. E. & Racine, J. S. (2012) Jackknife model averaging.
    \emph{Journal of Econometrics}, 979: 38–46
}
\seealso{
\code{\link{Weights}}, \code{\link{model.avg}}

Other model.weights: \code{\link{BGWeights}},
  \code{\link{bootWeights}}, \code{\link{cos2Weights}},
  \code{\link{stackingWeights}}
}
\keyword{models}

