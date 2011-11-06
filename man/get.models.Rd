\name{get.models}
\alias{get.models}
\encoding{utf-8}

\title{Get models}

\description{
    Generate a list of fitted model objects from a \code{model.selection} table.
}


\usage{
get.models(object, subset = delta <= 4, ...)
}

\arguments{
  \item{object}{ object returned by \code{\link{dredge}}. }
  \item{subset}{subset of models, an expression evaluated within the model
  selection table, see \link[=subset.model.selection]{subset method}. }
  \item{\dots}{additional parameters passed to \code{update}. For example, in
    \code{lme/lmer} one may want to use \code{method = "REML"} while using "ML"
    for model selection.}
}

\value{
    \code{\link{list}} of fitted model objects
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{dredge}}, \code{\link{model.avg}}
}

\examples{
# Mixed models:

require(nlme)
fm2 <- lme(distance ~ age + Sex, data = Orthodont,
    random = ~ 1 | Subject, method="ML")
ms2 <- dredge(fm2)

# Get top-most models, but fitted by REML:
(confset.d4 <- get.models(ms2, subset = delta < 4, method = "REML"))

}

\keyword{models}
