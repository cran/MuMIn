\name{AICc}
\alias{AICc}
\encoding{utf-8}
\title{Second-order Akaike Information Criterion}
\description{
Calculates second-order Akaike information criterion for one or several fitted
model objects (AIC for small samples). }

\usage{
AICc(object, ..., k = 2, REML = NULL)
}

\arguments{
  \item{object}{a fitted model object for which there exists a \code{logLik}
    method, or a \code{logLik} object}
  \item{\dots}{optionally more fitted model objects}
  \item{k}{the \dQuote{penalty} per parameter to be used; the default \code{k} =
    2 is the classical \code{\link[stats]{AIC}}}
   \item{REML}{ optional logical value, passed to the \code{logLik} method
    indicating whether the restricted log-likelihood or log-likelihood should be
    used. The default is to use the method used for model estimation.
   }
}

\value{
    If just one object is provided, returns a numeric value with the
    corresponding AICc; if more than one object are provided, returns a
    \code{data.frame} with rows corresponding to the objects and columns
    representing the number of parameters in the model (df) and AICc.
}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed. }


\author{Kamil Barto\enc{ń}{n}}

\seealso{
Akaike's An Information Criterion: \code{\link[stats]{AIC}}

\code{\link[AICcmodavg]{AICc}} in package \pkg{AICcmodavg},
\code{\link[glmulti]{aicc}} in package \pkg{glmulti}

}

\examples{

#Model-averaging mixed models

library(nlme)
data(Orthodont, package = "nlme")

# Fit model by REML
fm2 <- lme(distance ~ Sex*age, data = Orthodont,
    random = ~ 1|Subject / Sex, method = "REML")

# Model selection: ranking by AICc using ML
dd <- dredge(fm2, trace=TRUE, rank="AICc", REML=FALSE)

(attr(dd, "rank.call"))

# Get the models (fitted by REML, as in the global model)
gm <- get.models(dd, 1:4)

# Because the models originate from 'dredge(..., rank=AICc, REML=FALSE)',
# the default weights in 'model.avg' are ML based:
model.avg(gm, method = "NA")
# same result
#model.avg(gm, method = "NA", rank="AICc", rank.args = list(REML=FALSE))
# REML based weights
model.avg(gm, method = "NA", rank="AICc", rank.args = list(REML=TRUE))

}


\keyword{models}
