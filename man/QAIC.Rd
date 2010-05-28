\name{QAIC}
\alias{QAIC}
\encoding{utf-8}
\title{Quasi AIC}
\description{
Calculates \dQuote{quasi AIC} for one or several fitted model objects. This function is provided just as an example of custom rank function for use with \code{\link{model.avg}} and \code{\link{dredge}}
}

\usage{
QAIC(object, ..., chat)
}

\arguments{
  \item{object}{a fitted model object.}
  \item{\dots}{ optionally more fitted model objects.}
  \item{chat}{c - hat}
}

\value{
	If just one object is provided, returns a numeric value with the corresponding QAIC; 
	if more than one object are provided, returns a data.frame with rows corresponding to the objects.
}

\details{
\code{rank} is specified as a function or a symbol (e.g. a backquoted name) or a character string specifying a function.

Function \code{rank} must be able to accept model as a first argument and must always return a scalar.
}


\author{ Kamil Bartoń}

\examples{
budworm <- data.frame(ldose = rep(0:5, 2), numdead = c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16), sex = factor(rep(c("M", "F"), c(6, 6))))
budworm$SF <- cbind(budworm$numdead, 20 - budworm$numdead)
budworm.qlg <- glm(SF ~ sex*ldose, family = quasibinomial, data = budworm)

dd1 <- dredge(budworm.qlg, rank = "QAIC", chat = summary(budworm.qlg)$dispersion)
gm1 <- get.models(dd1, 1:4)

model.avg(gm1)

model.avg(gm1[[1]], gm1[[2]], rank = "QAIC", rank.args = list(chat = 1))

}

\keyword{models}