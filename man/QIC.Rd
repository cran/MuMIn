\name{QIC}
\alias{QIC}
\alias{quasiLik}
\encoding{utf-8}
\title{QIC and Quasi-Likelihood for GEE}
\description{
Calculate quasi-likelihood under the independence model criterion (QIC)
  for Generalized Estimating Equations.
}

\usage{
QIC(object, ..., typeR = FALSE)
quasiLik(object, ...)
}

\arguments{
  \item{object}{a fitted model object of class \code{gee}, \code{geepack} or
        \code{yags}.}
  \item{\dots}{ for QIC, optionally more fitted model objects. }  
  \item{typeR}{ whether to calculate QIC(R) rather than QIC(I). }  
}

\value{
    If just one object is provided, returns a numeric value with the
    corresponding QIC; if more than one object are provided, returns a
    \code{data.frame} with rows corresponding to the objects and columns
    representing QIC and QIC\eqn{_{u}}{u}.
}

\note{
This implementation is based partly on (revised) code from packages \pkg{yags}
and \pkg{ape}. The functions are still in experimental stage and should be used
with caution.
}

\references{
Pan W. (2001) Akaike's Information Criterion in Generalized Estimating Equations.
\emph{Biometrics} 57: 120-125

Hardin W. J., Hilbe, J. M. (2003) \emph{Generalized Estimating Equations}. CRC Press.
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
Methods exist for 
\code{\link[gee]{gee}} (package \pkg{gee}),
\code{\link[geepack]{geeglm}} (\pkg{geepack}), and
\code{\link[yags]{yags}} (\pkg{yags}).
\code{yags} and
\code{\link[ape]{compar.gee}} from package \pkg{ape} both provide QIC values.

}

\keyword{models}

\examples{

library(geepack)
data(ohio)

fm1 <- geeglm(resp ~ age * smoke, id = id, data = ohio,
    family = binomial, corstr = "exchangeable", scale.fix = TRUE)
fm2 <- update(fm1, corstr = "ar1")
fm3 <- update(fm1, corstr = "unstructured")

model.sel(fm1, fm2, fm3, rank = QIC)

\dontrun{
# same result:
    dredge(fm1, m.min = 3, rank = QIC, varying = list(
    corstr = list("exchangeable", "unstructured", "ar1")
    ))      
}

}
