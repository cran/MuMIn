\name{model.avg}
\alias{model.avg}
\alias{print.averaging}
\encoding{utf-8}

\title{Model averaging}
\description{
averages models according to an information criterion
}

\usage{
model.avg(m1, ..., beta = FALSE, method = c("0", "NA"), rank = NULL, rank.args = NULL, alpha = 0.05)
}

\arguments{
  \item{m1}{a fitted model object or a list of such objects}
  \item{\dots}{more fitted model objects}
  \item{beta}{logical, should standardized coefficients be returned rather than normal ones?}
  \item{method}{if set to "0" (default), terms missing in one model are assumed to be 0's, otherwise they are omitted from the weighted average}
  \item{rank}{custom rank function (information criterion) to use instead of \code{AICc}, e.g. \code{QAIC} or \code{BIC}, may be omitted if \code{m1} is a list returned by \code{dredge}. See \sQuote{Details}}
  \item{rank.args}{optional \code{list} of arguments for the \code{rank} function}
  \item{alpha}{significance level for calculatinq confidence intervals}

}

\value{
  model.avg returns a list with elements:
  \item{coefficients}{model coefficients}
  \item{variance}{variance of coefficients}
  \item{avg.model}{averaged model summary (\code{\link{data.frame}} with columns: coef - averaged coefficients, var - unconditional variance estimator, ase - adjusted standard error estimator, lci, uci - unconditional confidence intervals)}
  \item{relative.importance}{relative variable importances}
  \item{beta}{(logical) were standardized coefficients used?}

}

\details{
\code{rank} is found by a call to \code{match.fun} and typically is specified as a function or a symbol (e.g. a backquoted name) or a character string specifying a function to be searched for from the environment of the call to lapply.\cr
Function \code{rank} must be able to accept model as a first argument and must always return a scalar.
}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel inference: a practical information-theoretic approach}. 2nd ed. 
}

\author{Kamil Bartoń}

\seealso{
\code{\link{dredge}}, \code{\link{get.models}}.
\code{\link{QAIC}} has examples of using custom rank function.
}

\examples{
# Example from Burnham and Anderson (2002), page 100:
data(Cement)
lm1 <- lm(y ~ ., data = Cement)
dd <- dredge(lm1)
dd

#models with delta.aicc < 4
model.avg(get.models(dd, subset = delta < 4)) # get averaged coefficients

#or as a 95\% confidence set:
top.models <- get.models(dd, cumsum(weight) <= .95)

model.avg(top.models) # get averaged coefficients

#topmost model:
top.models[[1]]

## Not run:
# using BIC (Schwarz's Bayesian criterion) to rank the models
BIC <- function(x) AIC(x, k = log(length(residuals(x))))

model.avg(top.models, rank=BIC)



## End(Not run)

}

\keyword{models}
