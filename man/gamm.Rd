\name{gamm-wrapper}
\alias{gamm}
\alias{logLik.gamm}
\encoding{utf-8}

\title{Updateable \code{gamm}}
\description{
Enables updating of the model objects fitted by \code{gamm} and \code{gamm4}
from packages \pkg{mgcv} and \pkg{gamm4}.
}

\usage{
gamm(formula, random = NULL, ..., lme4 = inherits(random, "formula"))
}

\arguments{
  \item{formula, random, \dots}{ arguments passed to \code{gamm} or
	\code{gamm4}. }
  \item{lme4}{ logical, if \code{TRUE} \code{gamm4} is used rather than
	\code{gamm}. If \code{TRUE}, the \code{random} argument must be provided as
	a \code{formula}.
  }
}

\value{
	Depending on the value of the 'lme4' switch, either a \code{gamm} or
	\code{gamm4} fitted model object. The only difference from the original
	object is an addition of the \code{call} component.
}

\details{
This function is just a wrapper for \code{gamm} and \code{gamm4}. The only
purpose of it is to add a \code{call} component, that is not provided by
\code{gamm*} as such. It allows \code{update} on the returned object, so also
makes possible using it in model selection with \code{dredge}.

This is only a temporary workaround and it is likely be removed soon.
}

\note{
To assure \code{gamm} is called \emph{via} this wrapper, in case it is masked by
the original \code{gamm} from \pkg{mgcv} (when \pkg{MuMIn} was loaded after
\pkg{mgcv}), use \code{MuMIn::\link[MuMIn]{gamm}} syntax.


\code{logLik} method for \code{gamm} models returns log-likelihood value from the
\code{lme} component.
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link[mgcv]{gamm}} and \code{\link[gamm4]{gamm4}}
}


\keyword{models}
