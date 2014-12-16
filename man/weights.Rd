\name{Weights}
\alias{Weights}
\encoding{utf-8}

\title{Akaike weights}
\description{
Calculate or extract normalized model likelihoods (\sQuote{Akaike weights}).
}

\usage{
Weights(x)
}

\arguments{
	\item{x}{ a numeric vector of information criterion values such as AIC, or
	objects returned by functions like \code{AIC}. There are also methods for
	extracting Akaike weights from a \code{model.selection} or \code{averaging}
	objects. }

}

\value{
	a numeric vector of normalized likelihoods.
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{importance}}

\code{\link{weights}}, which extracts fitting weights from model objects.
}

\examples{
fm1 <- glm(Prop ~ dose, data = Beetle, family = binomial)
fm2 <- update(fm1, . ~ . + I(dose^2))
fm3 <- update(fm1, . ~ log(dose))
fm4 <- update(fm3, . ~ . + I(log(dose)^2))

round(Weights(AICc(fm1, fm2, fm3, fm4)), 3)
}

\keyword{models}
