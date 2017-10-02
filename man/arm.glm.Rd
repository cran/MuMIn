\name{arm.glm}
\alias{arm.glm}
\alias{armWeights}
\encoding{utf-8}
\title{Adaptive Regression by Mixing}
\description{
Combine all-subsets GLMs using the ARM algorithm.
Calculate ARM weights for a set of models.
}

\usage{
arm.glm(object, R = 250, weight.by = c("aic", "loglik"), trace = FALSE)

armWeights(object, ..., data, weight.by = c("aic", "loglik"), R = 1000)
}

\arguments{
  \item{object}{for \code{arm.glm}, a fitted \dQuote{global} \code{glm} object.
    For \code{armWeights}, a fitted \code{\link{glm}} object, or a
    \code{list} of such, or an \code{\link[=model.avg]{"averaging"}} object. }
  \item{\dots}{more fitted model objects. }
  \item{R}{number of permutations. }
  \item{weight.by}{indicates whether model weights should be calculated with AIC
    or log-likelihood. }
  \item{trace}{if \code{TRUE}, information is printed during the running of
    \code{arm.glm}. }
  \item{data}{a data frame in which to look for variables for use with 
    \link[=predict]{prediction}. If omitted, the fitted linear predictors are used.}
%% \item{seed}{optionally, the random seed. See \code{\link{set.seed}}.}
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
\code{arm.glm} returns an object of class \code{"averaging"} contaning only 
\dQuote{full} averaged coefficients. See \code{\link{model.avg}} for object 
description.

\code{armWeights} returns a numeric vector of model weights.
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

\code{\link{Weights}} for assigning new model weights to an \code{"averaging"}
object.

Other implementation of ARM algorithm: \code{arms} in (archived) package
\bold{MMIX}.

Other kinds of model weights: \code{\link{BGWeights}}, 
  \code{\link{bootWeights}}, 
  \code{\link{cos2Weights}}, \code{\link{jackknifeWeights}},
  \code{\link{stackingWeights}}.

}

\examples{
fm <- glm(y ~ X1 + X2 + X3 + X4, data = Cement)

summary(am1 <- arm.glm(fm, R = 15))

mst <- dredge(fm)

am2 <- model.avg(mst, fit = TRUE)

Weights(am2) <- armWeights(am2, data = Cement, R = 15)

# differences are due to small R:
coef(am1, full = TRUE)
coef(am2, full = TRUE)


}

\keyword{models}
