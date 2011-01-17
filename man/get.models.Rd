\name{get.models}
\alias{get.models}
\encoding{utf-8}

\title{Get models}

\description{
	Gets list of models from a \code{model.selection} object
}


\usage{
get.models(dd, subset = delta <= 4, ...)
}

\arguments{
  \item{dd}{ object returned by \code{\link{dredge}}}
  \item{subset}{subset of models}
  \item{\dots}{additional parameters passed to \code{update}, for example, in
	\code{lme/lmer} one may want to use \code{method = "REML"} while using "ML"
	for model selection}
}

\value{
	\code{\link{list}} of models.
}

\author{\enc{Kamil Bartoń}{Kamil Barton}}

\seealso{
\code{\link{dredge}}, \code{\link{model.avg}}
}

\examples{
# Mixed models:

require(nlme)
fm2 <- lme(distance ~ age + Sex, data = Orthodont,
	random = ~ 1 | Subject, method="ML")
dd2 <- dredge(fm2)

# Get top-most models, but fitted by REML:
(top.models.2 <- get.models(dd2, subset = delta < 4, method = "REML"))

}

\keyword{models}
