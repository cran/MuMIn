\name{MuMIn-package}
\alias{MuMIn-package}
\alias{MuMIn}
\docType{package}
\encoding{utf-8}
\title{Multi-model inference}
\description{
The package \code{MuMIn} contains functions to streamline model selection and
perform model averaging based on information criteria (AIC, AICc and alike).
}

\details{
User level functions include:
\describe{
	\item{\code{\link{dredge}}}{performs automated model selection with  
		subsets of the supplied \sQuote{global} model, and optional
		choices of other model properties (such as different link functions). 
		Model set may be generated either with \sQuote{all possible} 
		combinations, or tailored according to the conditions specified }
	\item{\code{\link{model.sel}}}{ creates a model selection table from 
		handpicked models }	
	\item{\code{\link{model.avg}}}{ calculates model averaged parameters, 
		with standard errors and confidence intervals }
	\item{\code{\link{AICc}}}{calculates second-order Akaike information 
		criterion.}
}

Model selection can be done according to any information criterion, such as 
\code{\link{AIC}}, \code{\link{AICc}}, \code{\link{BIC}}, \code{\link{QAIC}}, 
\code{\link{ICOMP}} or \link{Mallows' Cp}.

}

\author{ Kamil Barto\enc{ń}{n} }

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed.
}


\seealso{
\code{\link{AIC}}, \code{\link{step}} or \code{\link[MASS]{stepAIC}} for stepwise 
model selection by AIC.
}

\examples{

data(Cement)

fm1 <- lm(y ~ ., data = Cement)

ms1 <- dredge(fm1)
confset.d4 <- get.models(ms1, subset = delta < 4)
model.avg(confset.d4)

confset.95p <- get.models(ms1, cumsum(weight) <= .95)
avgmod.95p <- model.avg(confset.95p)
summary(avgmod.95p)
confint(avgmod.95p)

}

\keyword{package}
\keyword{models}
