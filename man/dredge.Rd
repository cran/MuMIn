\name{dredge}
\alias{dredge}
\alias{print.model.selection}

\encoding{utf-8}
\title{Automatized model selection}
\description{
Automatically generate models with combinations of the terms in the global
model, with optional rules for inclusion.
}

\usage{
dredge(global.model, beta = FALSE, evaluate = TRUE, rank = "AICc",
    fixed = NULL, m.max = NA, m.min = 0, subset, marg.ex = NULL,
    trace = FALSE, varying = NULL, ...)

\method{print}{model.selection}(x, abbrev.names = TRUE, ...)

}

\arguments{
	\item{global.model}{a fitted \sQuote{global} model object. See
	\sQuote{Details} for a list of supported types.
    }

	\item{beta}{logical, should standardized coefficients be returned?}

	\item{evaluate}{whether to evaluate and rank the models. If \code{FALSE}, a
	list of model \code{call}s is returned.}

	\item{rank}{optional custom rank function (information criterion) to be used
	instead \code{AICc}, e.g. \code{QAIC} or \code{BIC}. See \sQuote{Details}.}

	\item{fixed}{optional, either a single sided formula or a character vector
	giving names of terms to be included in all models.}

	\item{m.max, m.min}{optional, maximum and minimum number of terms in a single
	model (excluding the intercept), \code{m.max} defaults to the number of terms
	in \code{global.model}.}

	\item{subset}{ logical expression indicating models to keep in the resulting
		set. See \sQuote{Details}. }

	\item{marg.ex}{ a character vector specifying names of variables for which
	NOT to check for marginality restrictions when generating model formulas. If
	this argument is set to \code{TRUE}, all combinations of terms are used
	(i.e. no checking).}

	\item{trace}{
		if \code{TRUE}, all calls to the fitting function (i.e. updated
		\code{global.model} calls) are printed. }

	\item{varying}{ optionally, a named list describing the additional arguments
		to vary between the generated models. Names are the names of the
		arguments, and each item provides a list of choices. Complex items in
		the choice list (such as \code{family} objects) should be either named
		(uniquely) or quoted (unevaluated, e.g. using \code{\link{alist}}, see
		\code{\link{quote}}), otherwise it may produce rather unpleasant
		effects. See example in \code{\link{Beetle}. }
	}


  \item{x}{
    a \code{model.selection} object, returned by \code{dredge}.    }

  \item{abbrev.names}{ should variable names be abbreviated when printing?
  (useful with many variables). }

  \item{\dots}{ optional arguments for the \code{rank} function. Any can be
    an expression (of mode \code{call}), in which case any \code{x} within it
    will be substituted with a current model. }
}


\details{

\code{dregde} currently is known to work with \code{lm}, \code{glm},
\code{rlm}, \code{polr}, \code{multinom}, \code{gam}, \code{gls}, \code{lme},
\code{lmer}, \code{coxph}, \code{glmmML}, \code{sarlm}, \code{spautolm},
\code{gamm} and \code{gamm4} (the last two should be fitted \emph{via} the
wrapper \code{MuMIn::\link[MuMIn]{gamm}}).

Models are run one by one by repeated evaluation of modified calls to
\code{global.model} (in a similar fashion to that of \code{update}. This method,
while robust in that it can be applied to a variety of different model ofjects
is not very efficient, and may be time-intensive.

Note that the number of combinations grows exponentially with number of
predictor variables (\ifelse{latex}{\eqn{2^{N}}}{2^N}). Because there is
potentially a large number of models to evaluate, to avoid memory overflow the
fitted model objects are not stored. To get (a subset of) the models, use
\code{\link{get.models}} on the object returned by \code{dredge}.

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


The argument \code{subset} acts in a similar fashion to that in the function
\code{subset} for \code{data.frames}: the model terms can be referred to by name
as variables in the expression, with the difference that they are always logical
 (i.e. \code{TRUE} if a term exists in the model). The expression can contain
 any of the \code{global.model} terms (use \code{getAllTerms(global.model)} to
 list them). It can have a form of an unevaluated \code{call}, \code{expression}
 object, or a one sided \code{formula}. See \sQuote{Examples}.
Compound model terms (such as \sQuote{as-is} expressions within \code{I()} or
the smooths in \code{gam}) should be treated as non-syntactic names and
enclosed in back-ticks (see \link[base]{Quotes}). Mind the spacing, names
must match exactly the term names in model's formula. To simply keep certain
variables in all models, use of \code{fixed} is preferred.

Use of \code{na.action = na.omit} (\R's default) in \code{global.model} should
be avoided, as it results with sub-models fitted to different data sets, if
there are missing values. In versions >= 0.13.17 a warning is given in such a
case.
}


\value{
\code{dredge} returns an object of class \code{model.selection}, being a
\code{data.frame} with models' coefficients (or \code{TRUE}/\code{FALSE}
for factors), k, deviance/RSS, R-squared, AIC, AICc, delta and
 \emph{Akaike weight}. This depends on a type of model. Models are ordered
 according to the used information criterion (lowest on top), specified
 by \code{rank}.

The attribute \code{"calls"} is a list containing the model calls used (arranged
in the same order as the models).

}

\author{Kamil Barto\enc{ń}{n}}

\note{
Users should keep in mind the hazards that a \dQuote{thoughtless approach}
of evaluating all possible models poses. Although this procedure is in certain
cases useful and justified, it may result in selecting a spurious \dQuote{best}
model, due to model selection bias.

\emph{\dQuote{Let the computer find out} is a poor strategy and usually reflects
the fact that the researcher did not bother to think clearly about the problem
of interest and its scientific setting} (Burnham and Anderson, 2002).

}

\seealso{
\code{\link{get.models}}, \code{\link{model.avg}}.

There are \code{\link[=subset.model.selection]{subset}} and \code{plot} methods.

Possible alternatives: \code{\link[glmulti]{glmulti}} in package \pkg{glmulti}
and \code{\link[bestglm]{bestglm}} (\pkg{bestglm}),
or \code{\link[AICcmodavg]{aictab}} (\pkg{AICcmodavg}) and
\code{\link[bbmle]{ICtab}} (\pkg{bbmle}) for "hand-picked" model selection tables.

\code{\link[leaps]{regsubsets}} in package \pkg{leaps} also performs all-subsets
regression.


}


\examples{
# Example from Burnham and Anderson (2002), page 100:
data(Cement)
fm1 <- lm(y ~ ., data = Cement)
dd <- dredge(fm1)
subset(dd, delta < 4)

# Visualize the model selection table:
if(require(graphics))
plot(dd)


# Model average models with delta AICc < 4
model.avg(get.models(dd, subset = delta < 4))

#or as a 95\% confidence set:
confset.95p <- get.models(dd, cumsum(weight) <= .95)

model.avg(confset.95p) # get averaged coefficients

#'Best' model
summary(confset.95p[[1]])

\dontrun{
# Examples of using 'subset':
# exclude models containing both X1 and X2
dredge(fm1, subset = !(X1 & X2))
# keep only models containing X3
dredge(fm1, subset = ~ X3) # subset as a formula
dredge(fm1, subset = expression(X3)) # subset as expression object
# the same, but more effective:
dredge(fm1, fixed = "X3")

#Reduce the number of generated models, by including only those with
# up to 2 terms (and intercept)
dredge(fm1, m.max = 2)

}
}

\keyword{models}
