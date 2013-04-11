\name{get.models}
\alias{get.models}
\alias{pget.models}
\encoding{utf-8}

\title{Get models}

\description{
    Generate a list of fitted model objects from a \code{model.selection} table.
	\code{pget.models} can use paralell computation in a cluster to do that.
}


\usage{
get.models(object, subset, ...)
pget.models(object, cluster = NA, subset, ...)
}

\arguments{
  \item{object}{ object returned by \code{\link{dredge}}. }
  \item{subset}{subset of models, an expression evaluated within the model
  selection table, see the \link[=subset.model.selection]{subset method}.
	If it is a character vector, it is interpreted as names of rows to be
	selected. By default, \bold{all} model objects are fitted and returned 
	(see \sQuote{Note}). }
  \item{\dots}{additional arguments to update the models. For example, in
    \code{lme} one may want to use \code{method = "REML"} while using "ML"
    for model selection.}
	\item{cluster}{ a \code{cluster} object. See \code{\link{pdredge}} for
	details. }
}

\value{
    \code{\link{list}} of fitted model objects.
}

\note{
As of version 1.6.3, the default behaviour (if \code{subset} argument is
missing) is to return all the models, rather than a \sQuote{confidence set}
with \code{delta <= 4}.
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{dredge}}, \code{\link{model.avg}}
}

\examples{
# Mixed models:

require(nlme)
fm2 <- lme(distance ~ age + Sex, data = Orthodont,
    random = ~ 1 | Subject, method = "ML")
ms2 <- dredge(fm2)

# Get top-most models, but fitted by REML:
(confset.d4 <- get.models(ms2, subset = delta < 4, method = "REML"))

}

\keyword{models}
