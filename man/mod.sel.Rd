\name{model.sel}
\alias{mod.sel}
\alias{model.sel}
\alias{model.sel.default}
\alias{model.sel.model.selection}
\alias{model.sel}
\encoding{utf-8}

\title{model selection table}
\description{
Build a model selection table.
}

\usage{

model.sel(object, ...)

\method{model.sel}{default}(object, ..., rank = NULL, rank.args = NULL,
    beta = c("none", "sd", "partial.sd"), extra)
\method{model.sel}{model.selection}(object, rank = NULL, rank.args = NULL, fit = NA,
    ..., beta = c("none", "sd", "partial.sd"), extra)
}

\arguments{
	\item{object}{a fitted model object, a list of such objects, or a
	\code{"model.selection"} object.}

    \item{\dots}{more fitted model objects. }
	
    \item{rank}{optional, custom rank function (returning an information 
		criterion) to use instead of the default \code{AICc}, e.g. \code{QAIC} or \code{BIC}, may be omitted
        if \code{object} is a model list returned by \code{get.models}. }

    \item{rank.args}{optional \code{list} of arguments for the \code{rank}
        function. If one is an expression, an \code{x} within it is substituted
        with a current model. }
		
	\item{fit}{logical, stating whether the model objects should be re-fitted if
		they are not stored in the \code{"model.selection"} object. Set to
		\code{NA} to re-fit the models only if this is needed.
		See \sQuote{Details}. }
	
	\item{beta}{indicates whether and how the component models' coefficients 
        should be standardized. See the argument's description in 
        \code{\link{dredge}}. }
	
	\item{extra}{optional additional statistics to include in the result,
		provided as functions, function names or a list of such (best if named
		or quoted). See \code{\link{dredge}} for details. }

}

\value{
An object of class \code{"model.selection"} with columns containing useful
information about each model: the coefficients, \emph{df}, log-likelihood, the
value of the information criterion used, 
 \ifelse{latex}{\eqn{\Delta_{IC}}}{\enc{Δ}{Delta}_IC} and \sQuote{Akaike 
 weight}. If any arguments differ between the modelling function calls, the 
 result will include additional columns showing them (except for formulas and 
 some other arguments).
}

\details{
\code{model.sel} used with \code{"model.selection"} object will re-fit model
objects, unless they are stored in \code{object} (in attribute \code{"modelList"}),
 if argument \code{extra} is provided, or the requested \code{beta} is different
 than object's \code{"beta"} attribute, or the new \code{rank} function
 cannot be applied directly to \code{logLik} objects, or new \code{rank.args}
 are given (unless argument \code{fit = FALSE}).
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{dredge}}, \code{\link{AICc}}, \link[=MuMIn-models]{list of supported
models}.

Possible alternatives: \code{\link[bbmle]{ICtab}} (in package \pkg{bbmle}), or
 \code{\link[AICcmodavg]{aictab}} (\pkg{AICcmodavg}).
 
}

\examples{

Cement$X1 <- cut(Cement$X1, 3)
Cement$X2 <- cut(Cement$X2, 2)

fm1 <- glm(formula = y ~ X1 + X2 * X3, data = Cement)
fm2 <- update(fm1, . ~ . - X1 - X2)
fm3 <- update(fm1, . ~ . - X2 - X3)

## ranked with AICc by default
(msAICc <- model.sel(fm1, fm2, fm3))

## ranked with BIC
model.sel(fm1, fm2, fm3, rank = AIC, rank.args = alist(k = log(nobs(x))))
# or
# model.sel(msAICc, rank = AIC, rank.args = alist(k = log(nobs(x))))
# or
# update(msAICc, rank = AIC, rank.args = alist(k = log(nobs(x))))


}

\keyword{models}
