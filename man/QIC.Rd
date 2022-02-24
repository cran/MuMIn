\name{QIC}
\alias{QIC}
\alias{QICu}
\alias{quasiLik}
\encoding{utf-8}
\title{QIC and quasi-Likelihood for GEE}
\description{
Calculate quasi-likelihood under the independence model criterion (QIC)
  for Generalized Estimating Equations.
}

\usage{
QIC(object, ..., typeR = FALSE)
QICu(object, ..., typeR = FALSE)
quasiLik(object, ...)
}

\arguments{
  \item{object}{a fitted model object of class \code{"gee"}, \code{"geepack"}, 
    \code{"geem"}, \code{"wgee"},  or \code{"yags"}.}
  \item{\dots}{ for QIC and QIC\eqn{_{u}}{u}, optionally more fitted model objects. }  
  \item{typeR}{ logical, whether to calculate QIC(R). QIC(R) is
    based on quasi-likelihood of a working correlation \eqn{R} model. Defaults
    to \code{FALSE}, and QIC(I) based on independence model is returned.
        }
}

\value{
    If just one object is provided, returns a numeric value with the
    corresponding QIC; if more than one object are provided, returns a
    \code{data.frame} with rows corresponding to the objects and one column
    representing QIC or QIC\eqn{_{u}}{u}.
}

\note{
This implementation is based partly on (revised) code from packages \pkg{yags}
 (R-Forge) and \pkg{ape}.
}

\references{
Pan W. (2001) Akaike's Information Criterion in Generalized Estimating Equations.
\emph{Biometrics} 57: 120-125

Hardin J. W., Hilbe, J. M. (2003) \emph{Generalized Estimating Equations}.
Chapman & Hall/CRC
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
Methods exist for 
\code{\link[gee]{gee}} (package \pkg{gee}),
\code{\link[geepack]{geeglm}} (\pkg{geepack}), 
\code{\link[geeM]{geem}} (\pkg{geeM}), 
\code{\link[wgeesel]{wgee}} (\pkg{wgeesel}, the package's \code{QIC.gee} function is used), 
and
\code{yags} (\pkg{yags} on R-Forge).
There is also a \code{QIC} function in packages \pkg{MESS} and \pkg{geepack}, 
returning some extra information (such as CIC and QICc). \code{yags} and 
\code{\link[ape]{compar.gee}} from package \pkg{ape} both provide QIC values.
}

\keyword{models}

\examples{
\dontshow{ if(require(geepack)) \{ }

data(ohio)

fm1 <- geeglm(resp ~ age * smoke, id = id, data = ohio,
    family = binomial, corstr = "exchangeable", scale.fix = TRUE)
fm2 <- update(fm1, corstr = "ar1")
fm3 <- update(fm1, corstr = "unstructured")

# QIC function is also defined in 'geepack' but is returns a vector[6], so
# cannot be used as 'rank'. Either use `MuMIn::QIC` syntax or make a wrapper
# around `geepack::QIC`

QIC <- MuMIn::QIC
\dontrun{
QIC <- function(x) geepack::QIC(x)[1]
}

model.sel(fm1, fm2, fm3, rank = QIC)

#####
library(geepack)
library(MuMIn)

\dontrun{
# same result:
    dredge(fm1, m.lim = c(3, NA), rank = QIC, varying = list(
    corstr = list("exchangeable", "unstructured", "ar1")
    ))      
}
\dontshow{ \} }

}
