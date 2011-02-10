\name{dredge}
\alias{dredge}
\alias{print.model.selection}

\encoding{utf-8}
\title{Evaluate "all possible" models}
\description{
Automatically generate models with combinations of the terms in the global
model, with optional restrictions.
}

\usage{
dredge(global.model, beta = FALSE, eval = TRUE, rank = "AICc",
	fixed = NULL, m.max = NA, subset, marg.ex = NULL, trace = FALSE,
	...)

\method{print}{model.selection}(x, abbrev.names = TRUE, ...)

}

\arguments{
  \item{global.model}{a fitted \sQuote{global} model object. Currently, it can
	be a \code{lm}, \code{glm}, \code{rlm}, \code{multinom}, \code{gam},
	\code{gls}, \code{lme}, \code{lmer}, \code{sarlm} or \code{spautolm}, but
	also other types are likely to work (untested).
	}

  \item{beta}{logical, should standardized coefficients be returned?}

  \item{eval}{whether to evaluate and rank the models. If \code{FALSE}, a list
  of all possible model formulas is returned.}

  \item{rank}{optional custom rank function (information criterion) to be used
  instead \code{AICc}, e.g. \code{QAIC} or \code{BIC}. See \sQuote{Details}.}

  \item{fixed}{optional, either a single sided formula or a character vector
  giving names of terms to be included in all models.}

  \item{m.max}{optional, maximum number of terms to be included in single model,
  defaults to the number of terms in \code{global.model}.}

  \item{subset}{
	logical expression to put constraints for the set of models. Can contain any
	of the \code{global.model} terms (use \code{getAllTerms(global.model)} to
	list them). Complex expressions (e.g smooth functions in
	\code{\link[mgcv]{gam}} models) should be treated as non-syntactic names and
	enclosed in back-ticks (see \link[base]{Quotes}). Mind the spacing, names
	must match exactly the term names in model's formula.
	To simply keep certain variables in all models, use of \code{fixed} is
	preferred. }

  \item{marg.ex}{
	a character vector specifying names of variables for which NOT to check for
	marginality restrictions when generating model formulas. If this argument is
	set to \code{TRUE}, all model formulas are used (i.e. no checking). See
	\sQuote{Details}.}

  \item{trace}{
	if \code{TRUE}, all calls to the fitting function (i.e. updated
	\code{global.model} calls) are printed. }


  \item{x}{
	a \code{model.selection} object, returned by \code{dredge}.	}

  \item{abbrev.names}{
	Should variable names be abbreviated when printing? (useful with many
		variables). }

  \item{\dots}{optional arguments for the \code{rank} function. Any can be
	an expression (of mode \code{call}), in which case any \code{x} within it
	will be substituted with a current model. }
}


\details{

Models are run one by one by evaluating modified call of the \code{global.model}
\code{formula} argument (or \code{fixed} in \code{lme}). This method, while
robust in that it can be applied to a variety of different models is not very
efficient, and may be time consuming.

Because there is potentially a large number of models to evaluate, to avoid
memory overflow the fitted model objects are not stored. To get (a subset of)
the models, use \code{\link{get.models}} with the object returned by
\code{dredge} as an argument.

Handling interactions, \code{dredge} respects marginality constraints, so
\dQuote{all possible combinations} do not include models containing interactions
without their respective main effects. This behaviour can be altered by
\code{marg.ex} argument. It can be used to allow for simple nested designs. For
example, with global model of form \code{a / (x + z)}, use \code{marg.ex = "a"}
and \code{fixed = "a"}.

\code{rank} is found by a call to \code{match.fun} and may be specified as a
function or a symbol (e.g. a back-quoted name) or a character string specifying
a function to be searched for from the environment of the call to \code{dredge}.

Function \code{rank} must be able to accept model as a first argument and must
always return a scalar. Typical choice for \code{rank} would be "AIC", "QAIC" or
"BIC" (\pkg{stats} or \pkg{nlme}).


Use of \code{na.action = na.omit} (\R's default) in \code{global.model} should
be avoided, as it results with sub-models fitted to different data sets, if
there are missing values. In versions >= 0.13.17 a warning is given in such a
case.
}


\value{
\code{dredge} returns an object of class \code{model.selection}, being a
\code{\link{data.frame}} with models' coefficients (or \code{TRUE}/\code{FALSE}
for factors), k, deviance/RSS, R-squared, AIC, AICc, delta and weight. This
depends on a type of model. Models are ordered according to \code{\link{AICc}}
(lowest on top), or by \code{rank} function if specified.

The attribute "formulas" is a list containing model formulas (arranged in the
same order as the models).
}

\author{Kamil \enc{Bartoń}{Barton}}

\note{
Users should keep in mind the hazards that such a \dQuote{thoughtless approach}
of evaluating all possible models poses. Although this procedure is in certain
cases useful and justified, it may result in selecting a spurious \dQuote{best}
model, due to model selection bias.

\emph{\dQuote{Let the computer find out} is a poor strategy and usually reflects
the fact that the researcher did not bother to think clearly about the problem
of interest and its scientific setting} (Burnham and Anderson, 2002).
}

\seealso{
\code{\link{get.models}}, \code{\link{model.avg}}.
\code{\link{QAIC}} has examples of using custom rank function.

There is also \code{\link{subset.model.selection}} method.

Consider the alternatives: \code{\link[glmulti]{glmulti}} in package \pkg{glmulti}
and \code{\link[bestglm]{bestglm}} (\pkg{bestglm}),
or \code{\link[AICcmodavg]{aictab}} (\pkg{AICcmodavg}) and
\code{\link[bbmle]{ICtab}} (\pkg{bbmle}) for a "hand-picked" model selection tables.
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
# exclude models containing both X1 and X2
dredge(lm1, subset = !(X1 & X2))
# keep only models containing X3
dredge(lm1, subset = X3)
# the same, but more effective:
dredge(lm1, fixed = "X3")

#Reduce the number of generated models, by including only those with
# up to 2 terms (and intercept)
dredge(lm1, m.max = 2)

}
}

\keyword{models}
