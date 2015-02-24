\name{get.models}
\alias{get.models}
\alias{pget.models}
\encoding{utf-8}

\title{Retrieve models from selection table}

\description{
    Generate or extract a list of fitted model objects from a 
	\code{"model.selection"} table, optionally using parallel computation in a 
	cluster.
}


\usage{
get.models(object, subset, cluster = NA, ...)
}

\arguments{
  \item{object}{ object returned by \code{\link{dredge}}. }
  \item{subset}{subset of models, an expression evaluated within the model
	selection table (see \sQuote{Details}).}
	\item{cluster}{ optionally, a \code{"cluster"} object. If it is a valid 
		cluster, models are evaluated using parallel computation. }
	\item{\dots}{additional arguments to update the models. For example, in
    \code{lme} one may want to use \code{method = "REML"} while using \code{"ML"}
    for model selection.}
}

\value{
    \code{\link{list}} of fitted model objects.
}

\details{
The argument \code{subset} must be explicitely provided. This is to assure that 
a potentially long list of models is not fitted unintentionally. To evaluate all
 models, set \code{subset} to \code{NA} or \code{TRUE}. 

If \code{subset} is a character vector, it is interpreted as names of rows to be
selected. 

}

\note{
\code{pget.models} is still available, but is deprecated.
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{dredge}} and \code{\link{pdredge}}, \code{\link{model.avg}}

\code{\link{makeCluster}} in packages \pkg{parrallel} and \pkg{snow}

}

\examples{
# Mixed models:

\dontshow{ if(require(nlme)) \{ }
fm2 <- lme(distance ~ age + Sex, data = Orthodont,
    random = ~ 1 | Subject, method = "ML")
ms2 <- dredge(fm2)

# Get top-most models, but fitted by REML:
(confset.d4 <- get.models(ms2, subset = delta < 4, method = "REML"))

\dontrun{
# Get the top model:
get.models(ms2, subset = 1)[[1]]
}

\dontshow{  \} }
}

\keyword{models}
