\name{predict.averaging}
\alias{predict.averaging}
\encoding{utf-8}

\title{Predict Method for the Averaged Model}
\description{
Model-averaged predictions with optional standard errors. 
}

\usage{
\method{predict}{averaging}(object, newdata, se.fit = NULL, interval = NULL,
    type=c("link", "response"), ...)
}

\arguments{
    \item{object}{An object returned by \code{model.avg}.}
    \item{newdata}{ An optional data frame in which to look for variables with
        which to predict. If omitted, the fitted values are used.}
    \item{se.fit}{logical, indicates if standard errors should be returned.
        This has any effect only if the \code{predict} methods for each of the
        component models support it.}
    \item{interval}{Currently not used.}
    \item{type}{Predictions on response scale are only possible if all
        component models use the same \code{\link[stats]{family}}. See
        \code{\link[stats]{predict.glm}}}
    \item{\dots}{ arguments to be passed to respective \code{predict}
        method (e.g. \code{level} for \code{\link[nlme]{lme}} model). }

}

\value{

An object of class \code{averaging} with following elements:

 \item{summary}{a \code{data.frame} with deviance, AICc, Delta and weights
    for the component models.}
 \item{coefficients, variance}{matrices of component models' coefficients and
    their variances}
 \item{variable.codes}{names of the variables with numerical codes used in the
    summary}
 \item{avg.model}{the averaged model summary, (\code{data.frame} containing
    averaged coefficients, unconditional standard error, adjusted SE,
    and confidence intervals)}
 \item{importance}{the relative importance of variables}
 \item{beta}{(logical) were standardized coefficients used?}
 \item{term.names}{character vector giving names of all terms in the model}
 \item{residuals}{the residuals (response minus fitted values).}
 \item{x, formula}{the model matrix and formula analogical to those that would
    be used in a single model.}
 \item{method}{how the missing terms were handled ("NA" or "0").}
 \item{call}{the matched call.}

}

\details{

\code{predict.averaging} supports \code{method = "NA"} only for linear, fixed
    effect models. In other cases (e.g. nonlinear or mixed models), prediction
    is obtained using \dQuote{brute force}, i.e. by calling \code{predict} on
    each component model and weighted averaging the results, which is equivalent
    to assuming that missing coefficients equal zero (\code{method = "0"}).

Besides \code{predict} and \code{coef}, other generic methods such as
    \code{\link[stats]{formula}}, \code{\link[stats]{residuals}} and
    \code{\link[stats]{vcov}} are supported. \cr
    \code{logLik} method returns a list of \code{\link[stats]{logLik}} objects
    for the component models.
}

\note{

\code{predict.averaging} relies on availability of the \code{predict}
methods for the component model classes (except for \code{lm/glm}).

}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{model.avg}}
See \code{\link{par.avg}} for details of model-averaged parameter calculation.
}

\examples{

require(graphics)

# Example from Burnham and Anderson (2002), page 100:
data(Cement)
lm1 <- lm(y ~ ., data = Cement)
dd <- dredge(lm1)
top.models <- get.models(dd, subset=cumsum(weight) <= .95)
avgm <- model.avg(top.models)

# helper function
nseq <- function(x, len=length(x)) seq(min(x, na.rm=TRUE),
    max(x, na.rm=TRUE),    length=len)

# New predictors: X1 along the range of original data, other 
# variables held constant at their means
newdata <- as.data.frame(lapply(lapply(Cement[1:5], mean), rep, 25))
newdata$X1 <- nseq(Cement$X1, nrow(newdata))

# Predictions from each of the models in a set:
pred <- sapply(top.models, predict, newdata=newdata)
# Add predictions from the models averaged using two methods:
pred <- cbind(pred,
    averaged.0=predict(avgm, newdata),
    averaged.NA=predict(update(avgm, method="NA"), newdata))

matplot(x=newdata$X1, y=pred, type="l", lwd=c(rep(1,ncol(pred)-2), 2, 2),
    xlab="X1", ylab="y")

legend("topleft",
    legend=c(lapply(top.models, formula),
        paste("Averaged model (method=", c("0", "NA"), ")", sep="")),
    col=1:6, lty=1:5, lwd=c(rep(1,ncol(pred)-2), 2, 2), cex = .75)


\dontrun{
# Example with gam models (based on "example(gam)")
require(mgcv)
dat <- gamSim(1, n = 500, dist="poisson", scale=0.1)

gam1 <- gam(y ~ s(x0) + s(x1) + s(x2) +  s(x3) + (x1 + x2 + x3)^2,
    family = poisson, data = dat, method = "REML")

cat(dQuote(getAllTerms(gam1)), "\n")

# include only models with either smooth OR linear term (but not both)
# for each predictor variable:
dd <- dredge(gam1, subset=xor(`s(x1)`, x1) & xor(`s(x2)`, x2) & 
    xor(`s(x3)`, x3))
# ...this may take a while.

subset(dd, cumsum(weight) < .95)

top.models <- get.models(dd, cumsum(weight) <= .95)

newdata <- as.data.frame(lapply(lapply(dat, mean), rep, 50))
newdata$x1 <- nseq(dat$x1, nrow(newdata))
pred <- cbind(
    sapply(top.models, predict, newdata=newdata),
    averaged=predict(model.avg(top.models), newdata))

matplot(x=newdata$x1, y=pred, type="l", xlab="x1", ylab="y"
    lwd=c(rep(1, ncol(pred) - 2), 2, 2))
}


}

\keyword{models}
