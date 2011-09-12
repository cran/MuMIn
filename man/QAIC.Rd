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
    This implementation of QAIC uses model \code{\link{deviance}} rather than a likelihood itself. While this allows calculation also with models fitted using quasi-likelihood (where \code{logLik = NA}), the absolute values returned may differ from those obtained with the use of plain log-likelihood, since deviance is sometimes adjusted by a constant, so that the saturated model has deviance zero (see \code{\link{glm}}). \cr

    \code{\link{dredge}} will use \code{QAICc} instead of default \code{AICc} with \code{glm} with quasi* family.
}

\seealso{
\code{\link{AICc}} \cr
\code{\link[stats]{quasi}} family used for models with over-dispersion. \cr
\code{\link[stats]{AIC}} and \code{\link{BIC}} may also be used as a custom rank
function in \code{\link{dredge}} and \code{\link{model.avg}}.

\ifelse{html}{

\Sexpr[echo=TRUE,results=rd,stage=render]{
paste(
"\\\\sQuote{\\\\href{",
if(system.file(package="bbmle")!="") "../../bbmle/doc/" else "http://cran.r-project.org/web/packages/bbmle/vignettes/",
"quasi.pdf}", "{Dealing with quasi- models in R}}",
", a vignette in the \\\\pkg{bbmle} package.", sep="")
}
}{

\sQuote{\href{http://cran.r-project.org/web/packages/bbmle/vignettes/quasi.pdf}{Dealing with quasi- models in R}}, 
a vignette in the \pkg{bbmle} package.

}



}

%cat()


\author{Kamil Barto\enc{ń}{n}}

\examples{

# Based on "example(predict.glm)"
require(graphics)

budworm <- data.frame(
    ldose = rep(0:5, 2), 
    numdead = c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16), 
    sex = factor(rep(c("M", "F"), c(6, 6))))
budworm$SF = cbind(
    numdead = budworm$numdead, 
    numalive = 20 - budworm$numdead)

budworm.lg <- glm(SF ~ sex*ldose, data = budworm, family = quasibinomial)

dd <- dredge(budworm.lg, rank = "QAIC", 
    chat = summary(budworm.lg)$dispersion)
# Average all models
budworm.avg <- model.avg(get.models(dd, seq(nrow(dd))), method="NA")
#model.avg(mod[[1]], mod[[2]], rank = "QAIC", rank.args = list(chat = 1))

plot(c(1,32), c(0,1), type = "n", xlab = "dose",
     ylab = "prob", log = "x")
text(2^budworm$ldose, budworm$numdead/20, as.character(budworm$sex))
ld <- seq(0, 5, 0.1)

newdata <- data.frame(ldose=ld, sex=factor(rep("M", length(ld)),
    levels=levels(budworm$sex)))

# Predictions from global model / Males
pred.lg <- predict(budworm.lg, newdata, se.fit=TRUE, type="response")
matplot(2^ld, cbind(pred.lg$fit, pred.lg$fit - (2 * pred.lg$se.fit),
    pred.lg$fit + (2 * pred.lg$se.fit)), add=TRUE, type="l", col=1)

# Predictions from averaged model / Males
pred.avg <- predict(budworm.avg, newdata, se.fit=TRUE, type="response")
matplot(2^ld, cbind(pred.avg$fit, pred.avg$fit - (2 * pred.avg$se.fit),
    pred.avg$fit + (2 * pred.avg$se.fit)), add=TRUE, type="l", col=2)

newdata$sex[] <- "F"

# Predictions from global model / Females
pred.lg <- predict(budworm.lg, newdata, se.fit=TRUE, type="response")
matplot(2^ld, cbind(pred.lg$fit, pred.lg$fit - (2 * pred.lg$se.fit),
    pred.lg$fit + (2 * pred.lg$se.fit)), add=TRUE, type="l", col=1)

# Predictions from averaged model / Females
pred.avg <- predict(budworm.avg, newdata, se.fit=TRUE, type="response")
matplot(2^ld, cbind(pred.avg$fit, pred.avg$fit - (2 * pred.avg$se.fit),
    pred.avg$fit + (2 * pred.avg$se.fit)), add=TRUE, type="l", col=2)

legend("bottomright", legend=c("full", "averaged"), title="Model",
       col=1:2, lty=1)

}

\keyword{models}

