\name{AICc}
\alias{AICc}
\encoding{utf-8}
\title{Second-order Akaike Information Criterion}
\description{
Calculate Second-order Akaike Information Criterion for one or several fitted
model objects (AIC\eqn{_{c}}{c}, AIC for small samples). }

\usage{
AICc(object, ..., k = 2, REML = NULL)
}

\arguments{
  \item{object}{a fitted model object for which there exists a \code{logLik}
    method, or a \code{"logLik"} object. }
  \item{\dots}{optionally more fitted model objects. }
  \item{k}{the \sQuote{penalty} per parameter to be used; the default
  \code{k = 2} is the classical AIC. }
   \item{REML}{ optional logical value, passed to the \code{logLik} method
    indicating whether the restricted log-likelihood or log-likelihood should be
    used. The default is to use the method used for model estimation.
   }
}

\value{
If just one object is provided, returns a numeric value with the
corresponding AIC\eqn{_{c}}{c}; if more than one object are provided, returns a
\code{data.frame} with rows corresponding to the objects and columns
representing the number of parameters in the model (\emph{df}) and AIC\eqn{_{c}}{c}.
}

\note{
AIC\eqn{_{c}}{c} should be used instead AIC when sample size is small in
comparison to the number of estimated parameters (Burnham & Anderson 2002
recommend its use when \eqn{n / K < 40}).

}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed. New York,
Springer-Verlag.

Hurvich, C. M. and Tsai, C.-L. (1989) Regression and time series model selection
in small samples, \emph{Biometrika} 76: 297–307.
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
Akaike's An Information Criterion: \code{\link{AIC}}

Other implementations:
\code{\link[AICcmodavg]{AICc}} in package \pkg{AICcmodavg},
\code{\link[bbmle]{AICc}} in package \pkg{bbmle} and
\code{\link[glmulti]{aicc}} in package \pkg{glmulti}

}

\examples{
#Model-averaging mixed models
\dontshow{ if(require(nlme)) \{
oop <- }
options(na.action = "na.fail")

data(Orthodont, package = "nlme")

# Fit model by REML
fm2 <- lme(distance ~ Sex*age, data = Orthodont,
    random = ~ 1|Subject / Sex, method = "REML")

# Model selection: ranking by AICc using ML
ms2 <- dredge(fm2, trace = TRUE, rank = "AICc", REML = FALSE)

(attr(ms2, "rank.call"))

# Get the models (fitted by REML, as in the global model)
fmList <- get.models(ms2, 1:4)

# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:
summary(model.avg(fmList))

\dontrun{
# the same result:
model.avg(fmList, rank = "AICc", rank.args = list(REML = FALSE))
}
\dontshow{ \} }
}

\keyword{models}
