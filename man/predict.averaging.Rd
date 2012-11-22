\name{predict.averaging}
\alias{predict.averaging}
\encoding{utf-8}

\title{Predict method for the averaged model}
\description{
Model-averaged predictions with optional standard errors.
}

\usage{
\method{predict}{averaging}(object, newdata = NULL, se.fit = FALSE,
    interval = NULL, type = NA, backtransform = FALSE, full = TRUE, ...)
}

\arguments{
    \item{object}{An object returned by \code{model.avg}. }
    \item{newdata}{ An optional data frame in which to look for variables with
        which to predict. If omitted, the fitted values are used. }
    \item{se.fit}{logical, indicates if standard errors should be returned.
        This has any effect only if the \code{predict} methods for each of the
        component models support it. }
    \item{interval}{Currently not used. }
	\item{type}{The type of predictions to return (see documentation for
		\code{predict} appropriate for the class of used component models). If
		omitted, the default type is used. See \sQuote{Details}. }
    \item{backtransform}{If \code{TRUE}, the averaged predictions are
		back-transformed from link scale to response scale. This makes sense
		provided that all component models use the same \code{family}, and the
		prediction from each of the component models is calculated on the link 
		scale (as specified by \code{type}. For \code{glm}, use \code{type = 
		"link"}).
		See \sQuote{Details}.
		}
	\item{full}{ If \code{TRUE}, the full model averaged coefficients are used
		(only if \code{se.fit = FALSE} and the component objects are a result of
		\code{lm}). }
    \item{\dots}{Arguments to be passed to respective \code{predict}
        method (e.g. \code{level} for \code{\link[nlme]{lme}} model). }

}

\value{
If \code{se.fit = FALSE}, a vector of predictions, otherwise a list
with components: \code{fit} containing the predictions, and \code{se.fit} with
the estimated standard errors.
}

\details{
If all the component models are oridinary linear models, the prediction can be
made either with the full averaged coefficients (the argument \code{full =
TRUE} this is the default) or subset-averaged coefficients. Otherwise the
prediction is obtained by calling \code{predict} on each component model and
weighted averaging the results, which corresponds to the assumption that all
predictors are present in all models, but those not estimated are equal zero.
See \sQuote{Note} in \code{\link{model.avg}}. Predictions from component models
with standard errors are passed to \code{par.avg} and averaged in the same way
as the coefficients.

Predictions on the response scale from generalized models can be calculated by
averaging predictions of each model on the link scale, followed by inverse
transformation (this is achieved with \code{type = "link"} and
\code{backtransform = TRUE}). This is only possible if all component models use
the same family and link function. Alternatively, predictions from each model on
response scale may be averaged (with \code{type = "response"} and
\code{backtransform = FALSE}). Note that this leads to results differing from 
those calculated with the former method.
See also \code{\link[stats]{predict.glm}}.
}

\note{
This method relies on availability of the \code{predict} methods for the
component model classes (except when all component models are of class
\code{lm}).

The package \pkg{MuMIn} includes \code{predict} methods for \code{lme} and
\code{lmer} (\pkg{lme4}), both of which can calculate standard errors of the
predictions (with \code{se.fit = TRUE}).
The former enhances the \link[=predict.lme]{original function} from package
\pkg{nlme}, and with \code{se.fit = FALSE} it simply returns its result.
\pkg{MuMIn}'s version is always used in averaged model predictions (so it is
possible to predict with standard errors), but from within global environment it
will be found only if \pkg{MuMIn} is before \pkg{nlme} on the
\link[=search]{search list} (or directly extracted from namespace as
\code{MuMIn::predict.lme}).
\code{predict} method for \code{mer} models currently can only calculate values
on the outermost level (equivalent to \code{level = 0} in
\code{\link{predict.lme}}).
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
fm1 <- lm(y ~ X1 + X2 + X3 + X4, data = Cement)

ms1 <- dredge(fm1)
confset.95p <- get.models(ms1, subset = cumsum(weight) <= .95)
avgm <- model.avg(confset.95p)

nseq <- function(x, len = length(x)) seq(min(x, na.rm = TRUE),
    max(x, na.rm=TRUE), length = len)

# New predictors: X1 along the range of original data, other
# variables held constant at their means
newdata <- as.data.frame(lapply(lapply(Cement[1:4], mean), rep, 25))
newdata$X1 <- nseq(Cement$X1, nrow(newdata))

n <- length(confset.95p)

# Predictions from each of the models in a set, and with averaged coefficients
pred <- data.frame(
	model = sapply(confset.95p, predict, newdata = newdata),
	averaged.subset = predict(avgm, newdata, full = FALSE),
    averaged.full = predict(avgm, newdata, full = TRUE)
	)

opal <- palette(c(topo.colors(n), "black", "red", "orange"))
matplot(newdata$X1, pred, type = "l",
	lwd = c(rep(2,n),3,3), lty = 1,
    xlab = "X1", ylab = "y", col=1:7)

# For comparison, prediction obtained by averaging predictions of the component
# models
pred.se <- predict(avgm, newdata, se.fit = TRUE)
y <- pred.se$fit
ci <- pred.se$se.fit  * 2
matplot(newdata$X1, cbind(y, y - ci, y + ci), add = TRUE, type="l",
	lty = 2, col = n + 3, lwd = 3)

legend("topleft",
    legend=c(lapply(confset.95p, formula),
        paste(c("subset", "full"), "averaged"), "averaged predictions + CI"),
    lty = 1, lwd = c(rep(2,n),3,3,3),  cex = .75, col=1:8)

palette(opal)

}

\keyword{models}
