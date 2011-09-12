\name{MuMIn-package}
\alias{MuMIn-package}
\alias{MuMIn}
\docType{package}
\encoding{utf-8}
\title{
Multi-model inference
}
\description{
The package \code{MuMIn} contains functions for (automated) model selection and
model averaging based on information criteria (AIC alike).
}


\details{
User level functions include:
\describe{
    \item{\code{model.avg}}{does model averaging.}
	\item{\code{get.models}}{evaluates models from the table returned by dredge.}
	\item{\code{dredge}}{runs models with combinations of terms of the supplied
		\sQuote{global.model}.}
	\item{\code{AICc}}{calculates second-order Akaike information criterion for
		one or several fitted model objects.}
}
}

\author{ Kamil Barto\enc{ń}{n} \email{kamil.barton@go2.pl}}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed.
}


\seealso{
\code{\link[stats]{AIC}}, \code{\link[stats]{step}}

}
\examples{

fm1 <- lm(Fertility ~ . , data = swiss)

dd <- dredge(fm1)
top.models.1 <- get.models(dd, subset = delta < 4)
model.avg(top.models.1) # get averaged coefficients

top.models.2 <- get.models(dd, cumsum(weight) <= .95)
model.avg(top.models.2)

# Mixed models:
# modified example(lme)
data(Orthodont, package="nlme")
require(nlme)
fm2 <- lme(distance ~ age + Sex, data = Orthodont, 
    random = ~ 1 | Subject, method="ML")
dredge(fm2)

}

\keyword{package}
\keyword{models}
