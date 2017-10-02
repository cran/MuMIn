\name{Weights}
\alias{Weights}
\alias{Weights<-}
\encoding{utf-8}

\title{Akaike weights}
\description{
Calculate, extract or set normalized model likelihoods (\sQuote{Akaike weights}).
}

\usage{
Weights(x)
Weights(x) <- value
}

\arguments{
	\item{x}{a numeric vector of information criterion values such as AIC, or
	objects returned by functions like \code{AIC}. There are also methods for
	extracting \sQuote{Akaike weights} from \code{"model.selection"} or 
	\code{"averaging"} objects. }
	
	\item{value}{numeric, the new weights for the \code{"averaging"} object or
	\code{NULL} to reset the weights based on the original IC used.}

}

\details{
	The replacement function can assign new weights to an \code{"averaging"}
	object, affecting coefficient values and order of component models.
}

\value{
	For the extractor, a numeric vector of normalized likelihoods.
}


\note{
	On assigning new weights, the model order changes accordingly, so assigning
	the same weights again will cause incorrect re-calculation of averaged
	coefficients. To avoid that, either re-set model weights by assigning \code{NULL},
	or use ordered weights.
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{importance}}, \code{\link{weighted.mean}}

\code{\link{armWeights}}, 
\code{\link{bootWeights}}, \code{\link{BGWeights}}, \code{\link{cos2Weights}},
\code{\link{jackknifeWeights}} and \code{\link{stackingWeights}} can be used to
produce model weights.

\code{\link{weights}}, which extracts fitting weights from model objects.
}

\examples{

fm1 <- glm(Prop ~ dose, data = Beetle, family = binomial)
fm2 <- update(fm1, . ~ . + I(dose^2))
fm3 <- update(fm1, . ~ log(dose))
fm4 <- update(fm3, . ~ . + I(log(dose)^2))

round(Weights(AICc(fm1, fm2, fm3, fm4)), 3)


am <- model.avg(fm1, fm2, fm3, fm4, rank = AICc)

coef(am)

# Assign equal weights to all models:
Weights(am) <- rep(1, 4) # assigned weights are rescaled to sum to 1
Weights(am)
coef(am)

# Assign dummy weights:
wts <- c(2,1,4,3)
Weights(am) <- wts
coef(am)
# Component models are now sorted according to the new weights.
# The same weights assigned again produce incorrect results!
Weights(am) <- wts
coef(am) # wrong!
#
Weights(am) <- NULL # reset to original model weights
Weights(am) <- wts 
coef(am) # correct

}

\keyword{models}
