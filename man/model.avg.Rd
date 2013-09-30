\name{model.avg}
\alias{model.avg}
\alias{model.avg.default}
\alias{model.avg.model.selection}
\alias{print.averaging}
\encoding{utf-8}

\title{Model averaging}
\description{
Model averaging based on an information criterion.
}

\usage{

model.avg(object, ..., revised.var = TRUE)

\method{model.avg}{default}(object, ..., beta = FALSE, rank = NULL, rank.args = NULL,
    revised.var = TRUE, dispersion = NULL, ct.args = NULL)

\method{model.avg}{model.selection}(object, subset, fit = FALSE, ..., revised.var = TRUE)


}

\arguments{
	\item{object}{ a fitted model object or a list of such objects, or a
		\code{model.selection} object. See \sQuote{Details}. }

    \item{\dots}{ for default method, more fitted model objects. Otherwise,
		arguments that are passed to the default method. }

    \item{beta}{ logical, should standardized coefficients be returned? }

	\item{rank}{ optionally, a custom rank function (information criterion) to
		use instead of \code{AICc}, e.g. \code{BIC} or \code{QAIC}, may be
		omitted if \code{object} is a model list returned by \code{get.models}
		or a \code{model.selection} object. See \sQuote{Details}. }

    \item{rank.args}{ optional \code{list} of arguments for the \code{rank}
        function. If one is an expression, an \code{x} within it is substituted
        with a current model. }

    \item{revised.var}{ logical, indicating whether to use revised formula for
        standard errors. See \code{\link{par.avg}}. }

	\item{dispersion}{ the dispersion parameter for the family used. See
		\code{\link{summary.glm}}. This is used currently only with \code{glm},
		is silently ignored otherwise. }
	
	\item{ct.args}{ optional list of arguments to be passed to
	     \code{\link{coefTable}} (besides \code{dispersion}). }
		 
	\item{subset}{ see \code{\link[=subset.model.selection]{subset}} method for
		\code{model.selection} object. }

	\item{fit}{ if \code{TRUE}, the component models are fitted using
		\code{get.models}. See \sQuote{Details}.
	}
	


}

\value{

An object of class \code{averaging} is a list with components:

\item{summary}{ a \code{data.frame} with log-likelihood, IC, Delta(IC) and
	Akaike weights for the component models. }
\item{coef.shrinkage}{ a vector of full model-averaged coefficients, see
	\sQuote{Note}. }
\item{coefArray}{ an \code{array} of component models' coefficients, their
	standard errors, and degrees of freedom. }
\item{term.codes}{ names of the terms with numerical codes used in
	the \code{summary}. }
\item{avg.model}{ the model averaged parameters. A \code{data.frame} containing
    averaged coefficients, unconditional standard error, adjusted SE (if
	\emph{df}s are available) and z-values (coefficient and SE) and
	significance (assuming a normal error distribution). }
\item{importance}{ relative importance of the predictor variables (including
	interactions), calculated as a sum of the \emph{Akaike weights} over all of
	the models in which the parameter of interest appears. }
\item{term.names}{ character vector giving names of all terms in the model. }
\item{x, formula}{ the model matrix and formula corresponding to the one that
	would be used in a single model. \code{formula} contains only the averaged
	coefficients. }
%%\item{residuals}{ model averaged residuals (response minus fitted values). }
\item{call}{ the matched call. }

In addition, the object has following attributes:
\item{modelList}{ a list of component model objects. }
\item{beta}{ logical, were standardized coefficients used? }
\item{revised.var}{ if \code{TRUE}, the standard errors were calculated with the
	revised formula (See \code{\link{par.avg}}). }
}

\details{

\code{model.avg} may be used either with a list of models, or directly with a
	\code{model.selection} object (e.g. returned by \code{dredge}). In the
	latter case, the models from the model selection table are not evaluated
	unless the argument \code{fit} is set to \code{TRUE} or some additional
	arguments are present (such as \code{rank} or \code{dispersion}). This
	results in much faster calculation, but has certain drawbacks, because the
	fitted component model objects are not stored, and some methods (e.g.
	\code{predict}, \code{fitted}, \code{model.matrix} or \code{vcov}) would not
	be available with the returned object.
	Otherwise, \code{get.models} is called prior to averaging, and \dots are
	passed to it.

For a list of model types that are accepted see \link[=MuMIn-models]{list of
supported models}.


\code{rank} is found by a call to \code{\link[base]{match.fun}} and typically is
    specified as a function or a symbol (e.g. a back-quoted name) or a character
    string specifying a function to be searched for from the environment of the
    call to lapply. \code{rank} must be a function able to accept model as a
    first argument and must always return a scalar.

Several standard methods for fitted model objects exist for class
\code{averaging}, including
	\code{summary}, \code{\link[=predict.averaging]{predict}},
	\code{coef}, \code{confint}, \code{\link{formula}},
	%%\code{\link{residuals}},
	\code{\link{vcov}}.
The \code{coef} method a accepts argument \code{full}, if set to \code{TRUE} the
	full model-averaged coefficients are returned, rather than subset-averaged
	ones.
\code{logLik} returns a list of \code{\link[stats]{logLik}} objects
    for the component models.

}

\note{

The \sQuote{subset} (or \sQuote{conditional}) average only averages over the
	models where the parameter appears. An alternative, the \sQuote{full}
	average assumes that a variable is included in every model, but in some
	models the corresponding coefficient is set to zero. Unlike the
	\sQuote{subset average}, it does not have a tendency of biasing the value
	away from zero. It is, however, an unresolved issue how the variance of this
	estimate should be calculated, therefore the standard errors and confidence
	interval are returned only for the subset-averaged coefficients (as from
	version >= 1.5.0 argument \code{method} is no longer accepted).

Averaging models with different contrasts for the same factor would yield
	nonsense results, currently no checking for contrast consistency is done.

From version 1.0.1, \code{print} method provides only a concise output
	(similarly as for \code{lm}). To print a full summary of the results use
	\code{summary} function. Confidence intervals can be obtained with
	\code{\link{confint}}.

}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed. New York,
Springer-Verlag.
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{

See \code{\link{par.avg}} for more details of model averaged parameter
calculation.

\code{\link{dredge}}, \code{\link{get.models}} \cr
\code{\link{AICc}} has examples of averaging models fitted by REML.

\code{\link[AICcmodavg]{modavg}} in package \pkg{AICcmodavg}, and
\code{\link[glmulti]{coef.glmulti}} in package \pkg{glmulti} also perform model
averaging.
}

\examples{

# Example from Burnham and Anderson (2002), page 100:
data(Cement)
fm1 <- lm(y ~ ., data = Cement)
(ms1 <- dredge(fm1))

#models with delta.aicc < 4
summary(model.avg(ms1, subset = delta < 4))

#or as a 95\% confidence set:
avgmod.95p <- model.avg(ms1, cumsum(weight) <= .95)
confint(avgmod.95p)

\dontrun{
# The same result, but re-fitting the models via 'get.models'
confset.95p <- get.models(ms1, cumsum(weight) <= .95)
model.avg(confset.95p)

# Force re-fitting the component models
model.avg(ms1, cumsum(weight) <= .95, fit = TRUE)
# Models are also fitted if additional arguments are given
model.avg(ms1, cumsum(weight) <= .95, rank = "AIC")
}

\dontrun{
# using BIC (Schwarz's Bayesian criterion) to rank the models
BIC <- function(x) AIC(x, k = log(length(residuals(x))))
model.avg(confset.95p, rank = BIC)
# the same result, using AIC directly, with argument k
# 'x' in a quoted 'rank' argument is substituted with a model object
# (in this case it does not make much sense as the number of observations is
# common to all models)
model.avg(confset.95p, rank = AIC, rank.args = alist(k = log(length(residuals(x)))))
}

}

\keyword{models}
