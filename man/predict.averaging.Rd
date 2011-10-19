\name{predict.averaging}
\alias{predict.averaging}
\encoding{utf-8}

\title{Predict Method for the Averaged Model}
\description{
Model-averaged predictions with optional standard errors.
}

\usage{
\method{predict}{averaging}(object, newdata = NULL, se.fit = FALSE,
	interval = NULL, type = c("link", "response"), full = FALSE, ...)

}

\arguments{
    \item{object}{An object returned by \code{model.avg} }
    \item{newdata}{ An optional data frame in which to look for variables with
        which to predict. If omitted, the fitted values are used.}
    \item{se.fit}{logical, indicates if standard errors should be returned.
        This has any effect only if the \code{predict} methods for each of the
        component models support it.}
    \item{interval}{Currently not used }
    \item{type}{Predictions on response scale are only possible if all
        component models use the same \code{\link[stats]{family}}. See
        \code{\link[stats]{predict.glm}}. }
	\item{full}{ If \code{TRUE}, the full model averaged coefficients are used
		(only if \code{se.fit = FALSE} and the component objects are a result of
		\code{lm}). }
    \item{\dots}{Arguments to be passed to respective \code{predict}
        method (e.g. \code{level} for \code{\link[nlme]{lme}} model). }

}

\value{
If \code{se.fit = FALSE}, a vector of predictions. If \code{se = TRUE}, a list
with components: \code{fit} containing the predictions, and \code{se.fit} with
the estimated standard errors.
}

\details{
For an ordinary linear model, the predicted values can be averaged either using
the subset-averaged coefficients (if argument \code{full = FALSE}, the default),
or the full averaged coefficients.
In other cases the prediction is obtained by calling \code{predict} on each
component model and weighted averaging the results, which corresponds to
assuming that missing coefficients are averaged only over are averaged over a
subset of models wherein a variable appears. See \sQuote{Note} in
\code{\link{model.avg}}.

Predictions on the response scale from generalized models are calculated by
averaging predictions of each model on the link scale, followed by inverse
	transformation.

}

\note{
This method relies on availability of the \code{predict} methods for the
component model classes (except when all component models are of class \code{lm}).
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
fm1 <- lm(y ~ ., data = Cement)

ms1 <- dredge(fm1)
confset.95p <- get.models(ms1, subset=cumsum(weight) <= .95)
avgm <- model.avg(confset.95p)

nseq <- function(x, len=length(x)) seq(min(x, na.rm=TRUE),
    max(x, na.rm=TRUE), length=len)

# New predictors: X1 along the range of original data, other
# variables held constant at their means
newdata <- as.data.frame(lapply(lapply(Cement[1:5], mean), rep, 25))
newdata$X1 <- nseq(Cement$X1, nrow(newdata))

# Predictions from each of the models in a set:
pred <- sapply(confset.95p, predict, newdata=newdata)


pred <- cbind(pred,
    averaged.subset = predict(avgm, newdata),
    averaged.full = predict(avgm, newdata, full = TRUE))

matplot(x=newdata$X1, y=pred, type="l", lwd=c(rep(1,ncol(pred)-2), 2, 2),
    xlab="X1", ylab="y", col=c(3:7, 1,2), lty=c(1:5,1,1))

legend("topleft",
    legend=c(lapply(confset.95p, formula),
        paste(c("subset", "full"), "averaged")),
    col=c(3:7, 1,2), lty=c(1:5,1,1), lwd=c(rep(1,ncol(pred)-2), 2, 2), cex = .75)

}

\keyword{models}
