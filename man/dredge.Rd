\name{dredge}
\alias{dredge}
\alias{print.model.selection}

\encoding{utf-8}
\title{Evaluate "all possible" models}
\description{
Runs models with all possible combinations of the explanatory variables in the
supplied model. }

\usage{
dredge(global.model, beta = FALSE, eval = TRUE, rank = "AICc",
	fixed = NULL, m.max = NA, subset, ...)
}

\arguments{
  \item{global.model}{a fitted \sQuote{global} model object. Currently, it can
	be a \code{lm}, \code{glm}, \code{gam}, \code{lme}, \code{lmer},
	\code{sarlm} or \code{spautolm}, but also other types are likely to work
	(untested). }

  \item{beta}{logical should standardized coefficients be returned rather than
  normal ones?}

  \item{eval}{whether to evaluate and rank the models. If FALSE, a list of all
  possible model formulas is returned}

  \item{rank}{optional custom rank function (information criterion) to be used
  instead AICc, e.g. \code{QAIC} or \code{BIC}, See \sQuote{Details}}

  \item{fixed}{optional, either a single sided formula or a character vector
  giving names of terms to be included in all models}

  \item{m.max}{optional, maximum number of terms to be included in single model,
  defaults to the number of terms in \code{global.model}}

  \item{subset}{
	logical expression to put additional constraints for the set of models. Can
	contain any of the \code{global.model} terms. Run
	\code{getAllTerms(global.model)} to list all the terms. Complex expressions
	(e.g smooth functions in \code{\link[mgcv]{gam}} models) should be treated
	as non-syntactic names and enclosed in backticks (see \link[base]{Quotes}).
	Mind the spacing, names must match exactly the term names in model's
	formula.
  }

  \item{\dots}{optional arguments for the \code{rank} function. Any can be
	an expression (of mode \code{call}), in which case any \code{x} within it
	will be substituted with a current model. }
}


\details{

Models are run one by one by calling \code{\link{update}} with modified
\code{formula} argument. This method, while robust in that it can be applied to
a variety of different models, is not very efficient, so may be time (and
memory) consuming.

Handling interactions, \code{dredge} respects marginality constraints, so
\dQuote{all possible combinations} do not include models containing interactions
without their respective main effects.

\code{rank} is found by a call to \code{match.fun} and typically is specified as
a function or a symbol (e.g. a backquoted name) or a character string specifying
a function to be searched for from the environment of the call to lapply.

Function \code{rank} must be able to accept model as a first argument and must
always return a scalar.
}


\value{
\code{dredge} returns an object of class \code{model.selection}, being a a
\code{\link{data.frame}} with models' coefficients, k, deviance/RSS, R-squared,
AIC, AICc, delta and weight. This depends on a type of model. Models are ordered
according to \code{\link{AICc}} (lowest on top), or by \code{rank} function if
specified. The attribute "formulas" is a list containing model formulas.
}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed. }

\author{Kamil Bartoń}


\note{

Make sure there is no a \code{na.action} set to \code{na.omit} in
\code{global.model}. This can result with models fitted to different data sets,
if there are NA's present. \code{dredge} cannot handle nested model designs
(formulas such as \code{y ~ a/b}) properly.


}

\seealso{
\code{\link{get.models}}, \code{\link{model.avg}}.
\code{\link{QAIC}} has examples of using custom rank function.

There is also \code{\link{subset.model.selection}} method.
}


\examples{
# Example from Burnham and Anderson (2002), page 100:
data(Cement)
lm1 <- lm(y ~ ., data = Cement)
dd <- dredge(lm1)
subset(dd, delta < 4)

#models with delta.aicc < 4
model.avg(get.models(dd, subset = delta < 4)) # get averaged coefficients

#or as a 95\% confidence set:
top.models <- get.models(dd, cumsum(weight) <= .95)

model.avg(top.models) # get averaged coefficients

#topmost model:
top.models[[1]]

\dontrun{
# Examples of using 'subset':
# exclude models with with both X1 and X2
dredge(lm1, subset = !X1 | !X2)
# keep only models with X3
dredge(lm1, subset = X3)
# the same, but more effective:
dredge(lm1, fixed = ~ X3)
}
}

\keyword{models}
