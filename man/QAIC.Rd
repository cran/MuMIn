\name{QAIC}
\alias{QAIC}
\alias{QAICc}
\encoding{utf-8}
\title{Quasi AIC or AICc}
\description{
Calculate a modification of Akaike's Information Criterion for overdispersed
count data (or its version corrected for small sample,
\dQuote{quasi-AIC\eqn{_{c}}{c}}), for one or several fitted model objects.
}

\usage{
QAIC(object, ..., chat, k = 2, REML = NULL)
QAICc(object, ..., chat, k = 2, REML = NULL)
}

\arguments{
  \item{object}{a fitted model object.}
  \item{\dots}{ optionally, more fitted model objects.}
  \item{chat}{ \eqn{\hat{c}}{c-hat}, the variance inflation factor. }
  \item{k}{the \sQuote{penalty} per parameter. }
  \item{REML}{ optional logical value, passed to the \code{logLik} method
    indicating whether the restricted log-likelihood or log-likelihood should be
    used. The default is to use the method used for model estimation.
   }
}

\value{
    If only one object is provided, returns a numeric value with the
    corresponding QAIC or QAIC\eqn{_{c}}{c}; otherwise returns a
    \code{data.frame} with rows corresponding to the objects.
}

\note{
\eqn{\hat{c}}{c-hat} is the dispersion parameter estimated from the global
model, and can be calculated by dividing model's deviance by the number of
residual degrees of freedom.

In calculation of QAIC, the number of model parameters is increased by 1 to
account for estimating the overdispersion parameter. Without overdispersion,
\eqn{\hat{c} = 1}{c-hat = 1} and QAIC is equal to AIC.

Note that \code{glm} does not compute maximum-likelihood estimates in models
within the \emph{quasi-} family. In case it is justified, and with a proper
caution, a workaround could be used by \sQuote{borrowing} the \code{aic}
element from the corresponding \sQuote{non-quasi} family (see \sQuote{Example}).
}

\seealso{
\code{\link{AICc}}, \code{\link[stats]{quasi}} family used for models with
over-dispersion
}

\author{Kamil Barto\enc{ń}{n}}

\examples{
\dontshow{oop <- }
options(na.action = "na.fail")

# Based on "example(predict.glm)", with one number changed to create
# overdispersion
budworm <- data.frame(
    ldose = rep(0:5, 2), sex = factor(rep(c("M", "F"), c(6, 6))),
    numdead = c(10, 4, 9, 12, 18, 20, 0, 2, 6, 10, 12, 16))
budworm$SF = cbind(numdead = budworm$numdead,
    numalive = 20 - budworm$numdead)

budworm.lg <- glm(SF ~ sex*ldose, data = budworm, family = binomial)
(chat <- deviance(budworm.lg) / df.residual(budworm.lg))

dredge(budworm.lg, rank = "QAIC", chat = chat)
dredge(budworm.lg, rank = "AIC")

\dontrun{
# A 'hacked' constructor for quasibinomial family object, that allows for
# ML estimation
x.quasibinomial <- function(...) {
    res <- quasibinomial(...)
    res$aic <- binomial(...)$aic
    res
}
QAIC(update(budworm.lg, family = x.quasibinomial), chat = chat)
}
\dontshow{options(oop)}
}

\keyword{models}
