\name{importance}
\alias{importance}
\encoding{utf-8}

\title{Relative variable importance}
\description{
Sum of \sQuote{Akaike weights} over all models including the explanatory
variable.
}

\usage{
importance(x)
}

\arguments{
    \item{x}{either a list of fitted model objects, or a \code{"model.selection"}
    or \code{"averaging"} object. }

}

\value{
    a numeric vector of relative importance values, named as the predictor
    variables.
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{Weights}}

\code{\link{dredge}}, \code{\link{model.avg}}, \code{\link{model.sel}}
}

\examples{

# Generate some models
fm1 <- lm(y ~ ., data = Cement, na.action = na.fail)
ms1 <- dredge(fm1)

# Importance can be calculated/extracted from various objects:
importance(ms1)
\dontrun{
importance(subset(model.sel(ms1), delta <= 4))
importance(model.avg(ms1, subset = delta <= 4))
importance(subset(ms1, delta <= 4))
importance(get.models(ms1, delta <= 4))
}

# Re-evaluate the importances according to BIC
# note that re-ranking involves fitting the models again

# 'nobs' is not used here for backwards compatibility
lognobs <- log(length(resid(fm1)))

importance(subset(model.sel(ms1, rank = AIC, rank.args = list(k = lognobs)),
    cumsum(weight) <= .95))

# This gives a different result than previous command, because 'subset' is
# applied to the original selection table that is ranked with 'AICc'
importance(model.avg(ms1, rank = AIC, rank.args = list(k = lognobs),
    subset = cumsum(weight) <= .95))

}

\keyword{models}
