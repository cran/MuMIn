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
\item{x}{a numeric vector of any information criterion (such as 
     \AIC, \AICc, \QAIC, \BIC) values, or objects returned by functions like 
     \code{AIC}. There are also methods for extracting \sQuote{Akaike weights} from 
     \code{"model.selection"} or \code{"averaging"} objects. }
	
\item{value}{numeric, the new weights for the \code{"averaging"} object or
     \code{NULL} to reset the weights based on the original \IC used. 
      The assigned \code{value} need not sum to one, but if they are all zero, the 
      result will be invalid (\code{NaN}). }

}

\details{

\sQuote{Akaike weights}, \myequation{\omega_{i}}{X}{w_i}, of a model \var{i} can be 
interpreted as the probability that the model is the best (approximating) model
given the data and the set of all models considered. The weights are calculated as:

\mydequation{
\omega_i = \frac{\exp(\Delta_i/2)}{\sum_{r=1}^{R}\exp(\Delta_r/2)}
}{X}{
w_i = exp(Delta_i / 2) / (sum(exp(Delta_r / 2))
}

where \DeltaIC{i} is the \IC difference of the  \var{i}-th model relative to the
smallest \IC value in the set of \var{R} models.

The replacement version of \code{Weights} can assign new weights to an 
    \code{"averaging"} object, affecting coefficient values and the order of 
    component models. Upon assignment, the weights are normalised to sum to one.
    
}

\value{
	For the extractor, a numeric vector of normalized likelihoods.
}


\note{
Assigning new weights changes the model order accordingly, so reassigning
 weights to the same object must take this new order into account, otherwise 
 the averaged coefficients will be calculated incorrectly. To avoid this, either
 re-set the model weights by assigning \code{NULL}, or sort the new weights 
 using the (decreasing) order of the previously assigned weights.
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\lcode{sw}, \lcode{weighted.mean}

\lcode{armWeights}, 
\lcode{bootWeights}, \lcode{BGWeights}, \lcode{cos2Weights},
\lcode{jackknifeWeights} and \lcode{stackingWeights} can be used to
produce various kinds of model weights.

Not to be confused with \lcode{weights}, which extracts fitting weights 
from model objects.
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
