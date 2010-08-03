\name{model.avg}
\alias{model.avg}
\alias{print.averaging}
\alias{coef.averaging}
\alias{predict.averaging}
\encoding{utf-8}

\title{Model averaging}
\description{
Model averaging based on an information criterion.
}

\usage{
model.avg(m1, ..., beta = FALSE, method = c("0", "NA"), rank = NULL,
	rank.args = NULL, alpha = 0.05)

\method{coef}{averaging}(object, ...)

\method{predict}{averaging}(object, newdata, se.fit = NULL, interval = NULL,
	type=NULL, ...)
}


\arguments{
	\item{m1}{A fitted model object or a list of such objects. See
		\sQuote{Details}}

	\item{beta}{Logical, should standardized coefficients be returned rather
		than normal ones?}
	\item{method}{If set to \dQuote{0} (default), terms missing in one model are
		assumed to be 0's, otherwise they are omitted from the weighted
		average. See \sQuote{Details}.}
	\item{rank}{Custom rank function (information criterion) to use instead of
	\code{AICc}, e.g. \code{QAIC} or \code{BIC}, may be omitted if \code{m1} is
		a list returned by \code{dredge}. See \sQuote{Details}.}
	\item{rank.args}{Optional \code{list} of arguments for the \code{rank}
		function. If one is an expression, an \code{x} within it is substituted
		with a current model. }
	\item{alpha}{Significance level for calculatinq confidence intervals.}
	\item{object}{An object returned by \code{model.avg}.}
	\item{newdata}{ An optional data frame in which to look for variables with
		which to predict. If omitted, the fitted values are used.}
	\item{se.fit,interval}{Currently not used.}
	\item{type}{Ignored. Only predictions on the link scale are allowed.
		Warning is given if user tries something else here. }
	\item{\dots}{ for \code{model.avg} - more fitted model objects, for
		\code{predict} - arguments to be passed to respective \code{predict}
		method }

}

\value{
	An object of class \code{averaging}, being a list with elements:

	\item{summary}{Model table with deviance, AICc, Delta and weight.}
	\item{coefficients}{the model coefficients}
	\item{variance}{variance of coefficients}
	\item{avg.model}{averaged model summary (\code{\link{data.frame}} with
	  columns: coef - averaged coefficients, var - unconditional variance estimator,
	  ase - adjusted standard error estimator, lci, uci - unconditional confidence
	  intervals)}
	\item{relative.importance}{relative variable importances}
	\item{variable.codes}{Variable names with numerical codes used in the summary}
	\item{relative.importance}{Relative importance of variables}
	\item{weights}{}
	\item{beta}{(logical) were standardized coefficients used?}
	\item{model}{the model matrix, analogical to one that would be used in a
	  single model.}
	\item{residuals}{ the residuals (response minus fitted values).}

}

\details{

\code{model.avg} has been tested to work with \code{lm}, \code{glm}; \code{rlm}
(\pkg{MASS}); \code{gam} (\pkg{mgcv}); \code{lme}, \code{gls} (\pkg{nlme}),
\code{lmer} (\pkg{lme4}), as well as with \code{sarlm} and \code{spautolm}
(\pkg{spdep}). Other types are also likely to work, in particular when they
inherit from the supported classes.
See \sQuote{Details} section of the \link{Miscellaneous} page to see how to
provide support for other types of models.

\code{rank} is found by a call to \code{match.fun} and typically is specified as
a function or a symbol (e.g. a backquoted name) or a character string specifying
a function to be searched for from the environment of the call to lapply.

Function \code{rank} must be able to accept model as a first argument and must
always return a scalar.

\code{predict.averaging} supports method="NA" only for linear, fixed effect
models. In other cases (e.g. nonlinear or mixed models), prediction is obtained
using \dQuote{brute force}, i.e. by calling \code{predict} on each component
model and weighted averaging the results, which is equivalent to assuming that
missing coefficients equal zero (method="0").

Apart from \code{predict} and \code{coef}, other default methods, such as
 \code{formula} and \code{residuals} may be used.
}

\note{

\code{predict.averaging} relies on availability of the \code{predict} methods
for the component model classes (except for \code{(g)lm}).
}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed. }

\author{Kamil Bartoń}

\seealso{
\code{\link{dredge}}, \code{\link{get.models}}.
\code{\link{QAIC}} has examples of using custom rank function.
}

\examples{

require(graphics)

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

\dontrun{
# using BIC (Schwarz's Bayesian criterion) to rank the models
BIC <- function(x) AIC(x, k=log(length(residuals(x))))
mav <- model.avg(top.models, rank=BIC)
}


# Predicted values
nseq <- function(x, len=length(x)) seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE),
	length=len)

# New predictors: X1 along the range of original data, other variables held
#	constant at their means
newdata <- as.data.frame(lapply(lapply(Cement[1:5], mean), rep, 25))
newdata$X1 <- nseq(Cement$X1, nrow(newdata))

# Predictions from each of the models in a set:
pred <- sapply(top.models, predict, newdata=newdata)
# Add predictions from the models averaged using two methods:
pred <- cbind(pred,
	averaged.0=predict(model.avg(top.models, method="0"), newdata),
	averaged.NA=predict(model.avg(top.models, method="NA"), newdata)
	)

matplot(x=newdata$X1, y=pred, type="l", lwd=c(rep(1,ncol(pred)-2), 2, 2),
	xlab="X1", ylab="y")

legend("topleft",
	legend=c(lapply(top.models, formula),
		paste("Averaged model (method=", c("0", "NA"), ")", sep="")),
	col=1:6, lty=1:5, lwd=c(rep(1,ncol(pred)-2), 2, 2), cex = .75
	)


\dontrun{
# Example with gam models (based on "example(gam)")
require(mgcv)
dat <- gamSim(1, n = 500, dist="poisson", scale=0.1)

gam1 <- gam(y ~ s(x0) + s(x1) + s(x2) +  s(x3) + (x1+x2+x3)^2,
	family = poisson, data = dat, method = "REML")

cat(dQuote(getAllTerms(gam1)), "\n")

# include only models with smooth OR linear term (but not both) for each variable:
dd <- dredge(gam1, subset=(!`s(x1)` | !x1) & (!`s(x2)` | !x2) & (!`s(x3)` | !x3))
# ...this may take a while.

subset(dd, cumsum(weight) < .95)

top.models <- get.models(dd, cumsum(weight) <= .95)

newdata <- as.data.frame(lapply(lapply(dat, mean), rep, 50))
newdata$x1 <- nseq(dat$x1, nrow(newdata))
pred <- cbind(
	sapply(top.models, predict, newdata=newdata),
	averaged=predict(model.avg(top.models), newdata)
)

matplot(x=newdata$x1, y=pred, type="l", lwd=c(rep(1,ncol(pred)-2), 2, 2),
	xlab="x1", ylab="y")
}


}

\keyword{models}
