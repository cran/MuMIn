\name{model.avg}
\alias{model.avg}
\alias{print.averaging}
\encoding{utf-8}

\title{Model averaging}
\description{
Model averaging based on an information criterion.
}

\usage{

model.avg(object, ..., beta = FALSE, rank = NULL, rank.args = NULL,
    revised.var = TRUE)
}

\arguments{
    \item{object}{A fitted model object or a list of such objects. Alternatively
        an object of class \code{model.selection}. See \sQuote{Details}. }

    \item{\dots}{ more fitted model objects. }

    \item{beta}{Logical, should standardized coefficients be returned? }

    \item{rank}{Optional, custom rank function (information criterion) to use
        instead of \code{AICc}, e.g. \code{BIC} or \code{QAIC}, may be omitted
        if \code{object} is a model list returned by \code{get.models}
        or a \code{model.selection} object.
        See \sQuote{Details}. }

    \item{rank.args}{Optional \code{list} of arguments for the \code{rank}
        function. If one is an expression, an \code{x} within it is substituted
        with a current model. }

    \item{revised.var}{Logical, indicating whether to use revised formula for
        standard errors. See \code{\link{par.avg}}. }
}

\value{

An object of class \code{averaging} is a list with components:

\item{summary}{ a \code{data.frame} with log-likelihood, IC, Delta(IC) and 
	Akaike weights for the component models. }
\item{coefficients, se, dfs}{ matrices of component models' coefficients, their
	standard errors, and degrees of freedom. }
\item{coef.shrinkage}{ a vector of full model-averaged coefficients, see
	\sQuote{Note}. }
\item{variable.codes}{ names of the variables with numerical codes used in
    \code{summary}. }
\item{avg.model}{ the model averaged parameters. A \code{data.frame} containing
    averaged coefficients, unconditional standard error, adjusted SE (if
	\emph{df}s are available) and z-values (coefficient and SE) and
	significance (assuming a normal error distribution). }
\item{importance}{ relative importance of the predictor variables,
    calculated as a sum of the Akaike weights over all of the models in
    which the parameter of interest appears. }
\item{term.names}{ character vector giving names of all terms in the model. }
\item{x, formula}{ the model matrix and formula corresponding to the one that
	would be used in a single model. \code{formula} contains only the averaged
	coefficients. }
\item{residuals}{ model averaged residuals (response minus fitted values). }
\item{call}{ the matched call. }

In addition, the object has following attributes:
\item{mList}{ a list of component model objects. }
\item{beta}{ logical, were standardized coefficients used? }
\item{revised.var}{ if \code{TRUE}, the standard errors were calculated with the
	revised formula (See \code{\link{par.avg}}). }
}

\details{

\code{model.avg} has been tested to work with the fitted objects from the following
modelling functions:
	\code{lm}, \code{glm};
	\code{gam}, \code{gamm} (\pkg{mgcv});
	\code{gamm4} (\pkg{gamm4});
	\code{lme}, \code{gls} (\pkg{nlme});
	\code{lmer} (\pkg{lme4});
	\code{rlm}, \code{glm.nb}, \code{polr} (\pkg{MASS});
	\code{multinom} (\pkg{nnet});
	\code{sarlm}, \code{spautolm} (\pkg{spdep});
	\code{glmmML} (\pkg{glmmML});
	\code{coxph}, \code{survreg} (\pkg{survival});
	and several models within the class \code{unmarkedFit} (\pkg{unmarked}).
Other classes are also likely to be supported, in particular those
    inheriting from one of the above classes. See the vignette
	\sQuote{Extending \pkg{MuMIn}'s functionality} for a demonstration on how to
	provide support for other types of models.

\code{model.avg} may be used with a list of models, but also directly with a
	\code{model.selection} object returned by \code{dredge}. In the latter case,
	the models from the model selection table are evaluated (with a call to
	\code{get.models}) prior to averaging.
    A warning is given if the \code{subset} argument is not provided, and the
	default \code{delta <= 4} will be used.

\code{rank} is found by a call to \code{\link[base]{match.fun}} and typically is
    specified as a function or a symbol (e.g. a back-quoted name) or a character
    string specifying a function to be searched for from the environment of the
    call to lapply. \code{rank} must be a function able to accept model as a
    first argument and must always return a scalar.

Several standard methods for fitted model objects exist for class
\code{averaging}, including
	\code{summary}, \code{\link[=predict.averaging]{predict}},
	\code{coef}, \code{confint}, \code{\link{formula}},
	\code{\link{residuals}}, \code{\link{vcov}}.
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
	nonsense results.

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
summary(model.avg(get.models(ms1, subset = delta < 4))) # get averaged coefficients

#or as a 95\% confidence set:
confset.95p <- get.models(ms1, cumsum(weight) <= .95)

avgmod.95p <- model.avg(confset.95p) # get averaged coefficients
confint(avgmod.95p)


\dontrun{
# The same result
model.avg(ms1, cumsum(weight) <= .95)
}

\dontrun{
# using BIC (Schwarz's Bayesian criterion) to rank the models
BIC <- function(x) AIC(x, k=log(length(residuals(x))))
model.avg(confset.95p, rank=BIC)
# the same result, using AIC directly, with argument k
# 'x' in a quoted 'rank' argument is substituted with a model object
# (in this case it does not make much sense as the number of observations is
# common to all models)
model.avg(confset.95p, rank=AIC, rank.args=alist(k=log(length(residuals(x)))))
}

}

\keyword{models}
