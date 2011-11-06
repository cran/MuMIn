\name{MuMIn-package}
\alias{MuMIn-package}
\alias{MuMIn}
\docType{package}
\encoding{utf-8}
\title{Multi-model inference}
\description{
The package \pkg{MuMIn} contains functions to streamline information-theoretic 
model selection and carry out model averaging based on the information criteria.
}

\details{
The collection of functions includes:
\describe{
	\item{\code{\link{dredge}}}{performs automated model selection with
		subsets of the supplied \sQuote{global} model, and optional
		choices of other model properties (such as different link functions).
		The set of models may be generated either with \sQuote{all possible}
		combinations, or tailored according to the conditions specified. }
	\item{\code{\link{model.sel}}}{creates a model selection table from
		hand-picked models. }
	\item{\code{\link{model.avg}}}{ calculates model averaged parameters,
		with standard errors and confidence intervals. }
	\item{\code{\link{AICc}}}{calculates second-order Akaike information
		criterion. }
}

For a complete list of functions, use \code{library(help = "MuMIn")}. 

By default, AIC\eqn{_{c}}{c} is used to rank the models and to obtain model 
selection probabilities, though any other information criteria can be 
utilised. At least the following ones are currently implemented in \R: 
\code{\link{AIC}} and \code{\link{BIC}} in package \pkg{stats}, and
\code{\link{QAIC}}, \code{\link{QAICc}}, \code{\link{ICOMP}} and \link{Mallows' Cp} in \pkg{MuMIn}.

Most of \R's common modelling functions are supported, for a full list refer to the help pages for \code{dredge} and \code{model.avg}.
}

\author{ Kamil Barto\enc{ń}{n} }

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed. New York,
Springer-Verlag.
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
