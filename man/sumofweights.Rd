\name{sw}
\alias{importance}
\alias{sum.of.weights}
\alias{sw}
\encoding{utf-8}

\title{Per-variable sum of model weights}
\description{
Sum of model weights over all models including each explanatory variable.
}

\usage{
sw(x)
importance(x)
}

\arguments{
    \item{x}{either a list of fitted model objects, or a \code{"model.selection"}
    or \code{"averaging"} object. }

}

\value{
    a named numeric vector of so called relative importance values, for each 
    predictor variable.
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\lcode{Weights}

\lcode{dredge}, \lcode{model.avg}, \lcode{model.sel}
}

\examples{

# Generate some models
fm1 <- lm(y ~ ., data = Cement, na.action = na.fail)
ms1 <- dredge(fm1)

# Sum of weights can be calculated/extracted from various objects:
sw(ms1)
\dontrun{
sw(subset(model.sel(ms1), delta <= 4))
sw(model.avg(ms1, subset = delta <= 4))
sw(subset(ms1, delta <= 4))
sw(get.models(ms1, delta <= 4))
}

# Re-evaluate SW according to BIC
# note that re-ranking involves fitting the models again

# 'nobs' is not used here for backwards compatibility
lognobs <- log(length(resid(fm1)))

sw(subset(model.sel(ms1, rank = AIC, rank.args = list(k = lognobs)),
    cumsum(weight) <= .95))

# This gives a different result than previous command, because 'subset' is
# applied to the original selection table that is ranked with 'AICc'
sw(model.avg(ms1, rank = AIC, rank.args = list(k = lognobs),
    subset = cumsum(weight) <= .95))

}

\keyword{models}
