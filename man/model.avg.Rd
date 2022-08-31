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

\method{model.avg}{default}(object, ..., beta = c("none", "sd", "partial.sd"),
  rank = NULL, rank.args = NULL, revised.var = TRUE,
  dispersion = NULL, ct.args = NULL)

\method{model.avg}{model.selection}(object, subset, fit = FALSE, ..., revised.var = TRUE)


}

\arguments{
	\item{object}{ a fitted model object or a list of such objects, or a
		\code{"model.selection"} object. See \sQuote{Details}. }

    \item{\dots}{ for default method, more fitted model objects. Otherwise,
		arguments that are passed to the default method. }

    \item{beta}{ indicates whether and how the component models' coefficients 
        should be standardized. See the argument's description in 
        \code{\link{dredge}}. }

	\item{rank}{ optionally, a rank function (returning an information criterion) to
		use instead of \code{AICc}, e.g. \code{BIC} or \code{QAIC}, may be
		omitted if \code{object} is a model list returned by \code{get.models}
		or a \code{"model.selection"} object. See \sQuote{Details}. }

    \item{rank.args}{ optional \code{list} of arguments for the \code{rank}
        function. If one is an expression, an \code{x} within it is substituted
        with a current model. }

    \item{revised.var}{ logical, indicating whether to use the revised formula for
        standard errors. See \code{\link{par.avg}}. }

	\item{dispersion}{ the dispersion parameter for the family used. See
		\code{\link{summary.glm}}. This is used currently only with \code{glm},
		is silently ignored otherwise. }
	
	\item{ct.args}{ optional list of arguments to be passed to
	     \code{\link{coefTable}} (besides \code{dispersion}). }
		 
	\item{subset}{ see \code{\link[=subset.model.selection]{subset}} method for
		\code{"model.selection"} object. }

	\item{fit}{ if \code{TRUE}, the component models are fitted using
		\code{get.models}. See \sQuote{Details}.
	}
}

\value{

An object of class \code{"averaging"} is a list with components:

\item{msTable}{ a \code{data.frame} with log-likelihood, \emph{IC}, 
\ifelse{latex}{\eqn{\Delta_{IC}}}{\enc{Δ}{Delta}_IC} and
	\sQuote{Akaike weights} for the component models.
	Its attribute \code{"term.codes"} is a named vector with numerical 
    representation of the terms in the row names of \code{msTable}. 
}
\item{coefficients}{a \code{matrix} of model-averaged coefficients. 
    \dQuote{full}  coefficients in the first row, \dQuote{subset} coefficients in
    the second row. See \sQuote{Note} }
\item{coefArray}{ a 3-dimensional \code{array} of component models' coefficients,
	their standard errors and degrees of freedom. }
\item{sw}{ object of class \code{sw} containing per-model term 
    sum of model weights over all of the models in which the term appears.
}
\item{formula}{ a formula corresponding to the one that would be used in a 
	single model. The formula contains only the averaged (fixed) coefficients.
	}
\item{call}{ the matched call. }

The object has the following attributes:
\item{rank}{the rank function used. }
\item{modelList}{ optionally, a list of all component model objects. Only if
    the object was created with model objects (and not model selection table).}
\item{beta}{ Corresponds to the function argument. }
\item{nobs}{number of observations. }
\item{revised.var}{ Corresponds to the function argument. }

}

\details{

\code{model.avg} may be used either with a list of models or directly with a
	\code{model.selection} object (e.g. returned by \code{dredge}). In the
	latter case, the models from the model selection table are not evaluated
	unless the argument \code{fit} is set to \code{TRUE} or some additional
	arguments are present (such as \code{rank} or \code{dispersion}). This
	results in a much faster calculation, but has certain drawbacks, because the
	fitted component model objects are not stored, and some methods (e.g.
	\code{predict}, \code{fitted}, \code{model.matrix} or \code{vcov}) would not
	be available with the returned object. Otherwise, \code{get.models} is
	called prior to averaging, and \dots are passed to it.

For a list of model types that are accepted see \link[=MuMIn-models]{list of
	supported models}.

\code{rank} is found by a call to \code{\link[base]{match.fun}} and typically is
	specified as a function or a symbol or a character string specifying a 
	function to be searched for from the environment of the call to lapply. 
	\code{rank} must be a function able to accept model as a first argument and 
	must always return a numeric scalar.

Several standard methods for fitted model objects exist for class
	\code{averaging}, including \code{summary},
	\code{\link[=predict.averaging]{predict}}, \code{coef}, \code{confint},
	\code{\link{formula}}, and \code{\link{vcov}}.
	
\code{coef}, \code{vcov}, \code{confint} and \code{coefTable} accept argument
	\code{full} that if set to \code{TRUE}, the full model-averaged coefficients
	are returned, rather than subset-averaged ones (when \code{full = FALSE},
	being the default).

	
%	\code{
%	coef(object, full = FALSE, ...)
%	confint(object, parm, level = 0.95, full = FALSE, ...) 
%	vcov(object, full = FALSE, ...) 
%	coefTable(model, full = FALSE, ...)
%	}	
	
\code{logLik} returns a list of \code{\link[stats]{logLik}} objects
	for the component models.

}

\note{

The \sQuote{subset} (or \sQuote{conditional}) average only averages over the
	models where the parameter appears. An alternative, the \sQuote{full}
	average assumes that a variable is included in every model, but in some
	models the corresponding coefficient (and its respective variance) is set to
	zero.
	Unlike the \sQuote{subset average}, it does not have a tendency of
	biasing the value away from zero. The \sQuote{full} average is a type of
	shrinkage estimator, and for variables with a weak relationship to the
	response it is smaller than \sQuote{subset} estimators.
	
	Averaging models with different contrasts for the same factor would yield
	nonsense results. Currently, no checking for contrast consistency is done.
	
	\code{print} method provides a concise output (similarly as for \code{lm}). 
	To print more details use \code{summary} function, and \code{\link{confint}}
	to get confidence intervals.

}

\references{
Burnham, K. P. and Anderson, D. R. 2002 \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed. New York,
Springer-Verlag.

Lukacs, P. M., Burnham K. P. and Anderson, D. R. 2009 Model selection
bias and Freedman’s paradox. \emph{Annals of the Institute of 
Statistical Mathematics} \bold{62}, 117–125.

}

\author{Kamil Barto\enc{ń}{n}}

\seealso{

See \code{\link{par.avg}} for more details of model-averaged parameter
calculation.

\code{\link{dredge}}, \code{\link{get.models}} \cr
\code{\link{AICc}} has examples of averaging models fitted by REML.

\code{\link[AICcmodavg]{modavg}} in package \pkg{AICcmodavg}, and
\code{\link[glmulti]{coef.glmulti}} in package \pkg{glmulti} also perform model
averaging.
}

\examples{

# Example from Burnham and Anderson (2002), page 100:
fm1 <- lm(y ~ ., data = Cement, na.action = na.fail)
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
