\encoding{utf-8}
\name{bootWeights}
\alias{bootWeights}
\title{Bootstrap model weights}
\usage{
bootWeights(object, ..., R, rank = c("AICc", "AIC", "BIC"))
}
\arguments{
\item{object, \dots}{two or more fitted \code{\link{glm}} objects, or a
\code{list} of such, or an \code{\link[=model.avg]{"averaging"}} object.}

\item{R}{the number of replicates.}

\item{rank}{a character string, specifying the information criterion to use
for model ranking. Defaults to \code{\link{AICc}}.}
}
\value{
A numeric vector of model weights.
}
\description{
Compute model weights using bootstrap.
}
\details{
The models are fitted repeatedly to a resampled data set and ranked
using AIC-type criterion. The model weights represent the proportion of
replicates when a model has the lowest IC value.
}
\examples{
# To speed up the bootstrap, use 'x = TRUE' so that model matrix is included
#     in the returned object
fm <- glm(Prop ~ mortality + dose, family = binomial, data = Beetle, 
   na.action = na.fail, x = TRUE)

fml <- lapply(dredge(fm, eval = FALSE), eval)
am <- model.avg(fml)

Weights(am) <- bootWeights(am, data = Beetle, R = 25)

summary(am)

}
\references{
Dormann, C. et al. 2018 Model averaging in ecology: a review of Bayesian,
information-theoretic, and tactical approaches for predictive inference.
\emph{Ecological Monographs} \strong{88}, 485–504.
}
\seealso{
\code{\link{Weights}}, \code{\link{model.avg}}

Other model weights: 
\code{\link{BGWeights}()},
\code{\link{cos2Weights}()},
\code{\link{jackknifeWeights}()},
\code{\link{stackingWeights}()}
}
\author{
Kamil Barto\enc{ń}{n}, Carsten Dormann
}
\concept{model weights}
\keyword{models}
