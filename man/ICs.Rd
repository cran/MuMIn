\name{Information criteria}
\alias{IC}
\alias{Mallows' Cp}
\alias{Cp}
\alias{ICOMP}
\alias{CAICF}
\alias{DIC}
\encoding{utf-8}
\title{Various information criteria}
\description{
Calculate Mallows' \emph{Cp} and Bozdogan's ICOMP and CAIFC information criteria.

Extract or calculate Deviance Information Criterion from \code{MCMCglmm} and 
\code{merMod} object.
}

\usage{
Cp(object, ..., dispersion = NULL)
ICOMP(object, ..., REML = NULL)
CAICF(object, ..., REML = NULL)
DIC(object, ...)
}

\arguments{
  \item{object}{a fitted model object (in case of ICOMP and CAICF, \code{logLik}
	and \code{vcov} methods must exist for the object). For \code{DIC}, an 
	object of class \code{"MCMCglmm"} or \code{"merMod"}.
	}
  \item{\dots}{optionally more fitted model objects. }
  \item{dispersion}{the dispersion parameter. If \code{NULL}, it is inferred
	from object. }
  \item{REML}{ optional logical value, passed to the \code{logLik} method
    indicating whether the restricted log-likelihood or log-likelihood should be
    used. The default is to use the method used for model estimation. }
}

\value{
	If just one object is provided, the functions return a numeric value with
	the corresponding IC; otherwise a \code{data.frame} with rows corresponding
	to the objects is returned.
}

\details{

Mallows' \emph{Cp} statistic is the residual deviance plus twice the estimate of
\eqn{\sigma^{2}}{sigma^2} times the residual degrees of freedom. It is closely
related to AIC (and a multiple of it if the dispersion is known).

ICOMP (I for informational and COMP for complexity) penalizes the covariance
complexity of the model, rather than the number of parameters directly.

% \deqn{\mbox{ICOMP}=-2\log\mathit{Lik}(\hat{\theta})+k(\log(n)+2)+2C(\hat{\Sigma}_{model})}

%\deqn{\mbox{ICOMP}=-2\log L+k\,\log(\frac{tr(\Sigma)}{k})-\log\mid\Sigma\mid}
%
%where log\emph{L} is the log-likelihood, \eqn{\Sigma} is the variance-covariance
%matrix of parameters, \emph{k} is the number of parameters, and \emph{tr} is the
%trace of the matrix.

CAICF (C is for \sQuote{consistent} and F denotes the use of the Fisher
 information matrix) includes with penalty the natural logarithm of the
 determinant of the estimated Fisher information matrix.

% deqn{\mbox{CAICF}=-2\log\mathit{Lik}(\hat{\theta})+k(\log(n)+2)+
% \log|\mathit{I}(\hat{\theta})|}


}

% Bozdogan, H (1990) On the information-based measure of covariance complexity
% and its application to the evaluation of multivariate linear models. \emph{Comm.
% Stat. Theory and Methods} 19: 221-278
% Bozdogan, H. (2000). Akaike's Information Criteria and recent developments in
% Information Complexity. \emph{J. Math. Psych.} 44: 62-91

\references{
Mallows, C. L. 1973 Some comments on \emph{Cp}. \emph{Technometrics} \bold{15},
661–675.

Bozdogan, H. and Haughton, D. M. A. (1998) Information complexity criteria for
regression models. \emph{Comp. Stat. & Data Analysis} \bold{28}, 51--76.

Anderson, D. R. and Burnham, K. P. 1999 Understanding information criteria for selection among 
capture-recapture or ring recovery models. \emph{Bird Study} \bold{46}, 14--21.

Spiegelhalter, D. J., Best, N. G., Carlin, B. R., van der Linde, A. 2002 Bayesian measures of model 
complexity and fit. \emph{Journal of the Royal Statistical Society Series B-Statistical Methodology}
 \bold{64}, 583--616.

}

\seealso{
\code{\link{AIC}} and \code{\link{BIC}} in \pkg{stats}, \code{\link{AICc}}.
\code{\link{QIC}} for GEE model selection.
\code{\link[arm]{extractDIC}} in package \pkg{arm}, on which the (non-visible) 
method \code{extractDIC.merMod} used by \code{DIC} is based.
}

% \examples{}


\keyword{models}
