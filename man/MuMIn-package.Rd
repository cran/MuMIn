\name{MuMIn-package}
\alias{MuMIn-package}
\alias{MuMIn}
\docType{package}
\encoding{utf-8}
\title{
Multi model inference
}
\description{
The \code{MuMIn} contains functions for model selection and model averaging, based on information criteria (AIC alike).
}
\details{

User level functions include:\cr
\code{model.avg} - does model averaging.\cr
\code{get.models} - evaluates models from the table returned by dredge.\cr
\code{dredge} - runs models with combinations of terms of the supplied \sQuote{global} model.\cr
\code{AICc} - calculates second-order Akaike information criterion for one or several fitted model objects.\cr

}

\author{ Kamil Bartoń \email{kbarton@zbs.bialowieza.pl}}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel inference: a practical information-theoretic approach}. 2nd ed. 
}
\keyword{package}
\keyword{models}


\seealso{
\code{\link[stats]{AIC}}, \code{\link[stats]{step}}

}
\examples{

fm1 <- lm(Fertility ~ . , data = swiss)

dd <- dredge(fm1)
dd
#list of models with delta.aicc < 4
top.models.1 <- get.models(dd, subset = delta < 4)
model.avg(top.models.1) # get averaged coefficients

#or as a 95\% confidence set:
top.models.2 <- get.models(dd, cumsum(weight) <= .95)

model.avg(top.models.2) # get averaged coefficients

# Mixed models:

data(Orthodont, package="nlme")

require(nlme)
fm2 <- lme(distance ~ age + Sex, data = Orthodont, random = ~ 1 | Subject, method="ML")
dredge(fm2)


require(lme4)
fm3 <- lmer(distance ~ age + Sex + (1 | Subject), data = Orthodont, REML=FALSE)
dd3 <- dredge(fm3)

# Get top-most models, but fitted by REML:
(top.models.3 <- get.models(dd3, subset = delta < 4, REML=TRUE))
# use: method = "REML" for older versions of lme4

}