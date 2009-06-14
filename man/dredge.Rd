\name{dredge}
\alias{dredge}
\alias{print.model.selection}

\encoding{utf-8}
\title{data dredging}
\description{
runs models with all possible combinations of the explanatory variables in the supplied model.
}
\usage{
dredge(global.model, beta = FALSE, eval = TRUE, rank = "AICc",
		 fixed = NULL, m.max = NA, ...)
}

\arguments{
  \item{global.model}{a fitted \sQuote{global} model object. Currently, it can be a \code{lm}, \code{glm}, \code{gam}, \code{lme}, \code{lmer}, \code{sarlm} or \code{spautolm}.}
  \item{beta}{logical should standardized coefficients be returned rather than normal ones?}
  \item{eval}{whether to evaluate and rank the models. If set to FALSE only a list of all possible model formulas is returned}
  \item{rank}{custom rank function (information criterion) to use instead AICc, e.g. \code{QAIC} or \code{BIC}, See \sQuote{Details}}
  \item{fixed}{optional, either a single sided formula or a character vector giving names of terms to be included in all models}
  \item{m.max}{maximum number of terms to be included in single model, defaults to the number of terms in \code{global.model}}
  \item{\dots}{optional arguments for the \code{rank} function}
}


\details{
\code{rank} is found by a call to \code{match.fun} and typically is specified as a function or a symbol (e.g. a backquoted name) or a character string specifying a function to be searched for from the environment of the call to lapply. \cr

Function \code{rank} must be able to accept model as a first argument and must always return a scalar.
}


\value{
\code{dredge} returns a \code{\link{data.frame}} of class \code{model.selection} with models' coefficients, k, deviance/RSS, r-squared, AIC, AICc, etc. This depends on a type of model. Models are ordered according to \code{\link{AICc}} (lowest on top), or by \code{rank} function if specified. It has also a \code{formulas} attribute - a list containing all the model formulas.
  ...
}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel inference: a practical information-theoretic approach}. 2nd ed.
}

\author{Kamil Bartoń}


\note{
Use the \code{lmer} from \pkg{lme4} library rather than the old \code{lmer} from \pkg{Matrix} package.
Complex expressions (like \code{log(y)}) cannot be used as the response variable, the response must be specified as simple variable name. They may be used on the right side of the formula, though.
Make sure there is no a \code{na.action} set to \code{na.omit} in \code{global.model}. This can result with models fitted to different data sets, if there are NA's present.
\code{dredge} cannot handle nested model designs (formulas such as \code{y ~ a/b}).
}

\seealso{
\code{\link{get.models}}, \code{\link{model.avg}}.
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

}


\keyword{models}
