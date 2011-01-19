\name{QAIC}
\alias{QAIC}
\alias{QAICc}
\encoding{utf-8}
\title{Quasi AIC(c))}
\description{
Calculates \dQuote{quasi AIC} (or \dQuote{quasi AICc}) for one or several fitted model objects. This
function is provided mainly as an example of custom rank function for use with
\code{\link{model.avg}} and \code{\link{dredge}} }

\usage{
QAIC(object, ..., chat)
QAICc(object, ..., chat)
}

\arguments{
  \item{object}{a fitted model object.}
  \item{\dots}{ optionally more fitted model objects.}
  \item{chat}{c - hat}
}

\value{
	If just one object is provided, returns a numeric value with the
	corresponding QAIC; if more than one object are provided, returns a
	data.frame with rows corresponding to the objects.
}

\note{
This implementation of \code{QAIC}/\code{QAICc} may not be correct for mixed models. \cr
\code{\link{dredge}} will use \code{QAICc} instead of default \code{AICc} with \code{glm} with quasi*
family.
}

\seealso{
\code{\link{AICc}}

\code{\link[stats]{quasi}} family used for models with over-dispersion. \cr
\code{\link[stats]{AIC}}, \code{\link[nlme]{BIC}} in \pkg{nlme} (or \code{\link[stats4]{BIC}} in \pkg{stats4}) may
also be used as a custom rank function in \code{\link{dredge}} and \code{\link{model.avg}}.

}


\author{Kamil \enc{Bartoń}{Barton}}

\examples{

# Based on "example(predict.glm)"
require(graphics)

budworm <- data.frame(ldose = rep(0:5, 2), numdead = c(1, 4, 9, 13, 18, 20, 0,
	2, 6, 10, 12, 16), sex = factor(rep(c("M", "F"), c(6, 6))))
budworm$SF = cbind(numdead = budworm$numdead, numalive = 20 - budworm$numdead)

budworm.lg <- glm(SF ~ sex*ldose, data = budworm, family = quasibinomial)

plot(c(1,32), c(0,1), type = "n", xlab = "dose",
     ylab = "prob", log = "x")
text(2^budworm$ldose, budworm$numdead/20, as.character(budworm$sex))
ld <- seq(0, 5, 0.1)
lines(2^ld, predict(budworm.lg, data.frame(ldose=ld,
   sex=factor(rep("M", length(ld)), levels=levels(budworm$sex))),
   type = "response"))
lines(2^ld, predict(budworm.lg, data.frame(ldose=ld,
   sex=factor(rep("F", length(ld)), levels=levels(budworm$sex))),
   type = "response"))


dd <- dredge(budworm.lg, rank = "QAIC",
	chat = summary(budworm.lg)$dispersion)
mod <- get.models(dd, seq(nrow(dd)))
budworm.avg <- model.avg(mod)
model.avg(mod[[1]], mod[[2]], rank = "QAIC", rank.args = list(chat = 1))

linkinv <- quasibinomial()$linkinv
lines(2^ld, linkinv(predict(budworm.avg, data.frame(ldose=ld,
   sex=factor(rep("M", length(ld)), levels=levels(budworm$sex))))), col=2)
lines(2^ld, linkinv(predict(budworm.avg, data.frame(ldose=ld,
   sex=factor(rep("F", length(ld)), levels=levels(budworm$sex))))), col=2)
legend("bottomright", legend=c("full", "averaged"), title="Model",
	   col=1:2, lty=1)

}

\keyword{models}
