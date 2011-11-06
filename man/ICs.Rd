\name{Information criteria}
\alias{IC}
\alias{Mallows' Cp}
\alias{Cp}
\alias{ICOMP}
\encoding{utf-8}
\title{Various information criteria}
\description{
Calculate Mallows' \emph{Cp} and Bozdogan's ICOMP information criterion.
}

\usage{
Cp(object, dispersion = NULL)
ICOMP(object, ..., REML = NULL)
}

\arguments{
  \item{object}{a fitted model object (in case of ICOMP, a \code{logLik} method
	must exist for the object). }
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

ICOMP  (I for informational and COMP for complexity) penalizes the covariance
complexity of the model, rather than the number of parameters directly.

%\deqn{\mbox{ICOMP}=-2\log L+k\,\log(\frac{tr(\Sigma)}{k})-\log\mid\Sigma\mid}
%
%where log\emph{L} is the log-likelihood, \eqn{\Sigma} is the variance-covariance
%matrix of parameters, \emph{k} is the number of parameters, and \emph{tr} is the
%trace of the matrix.

}

% Bozdogan, H (1990) On the information-based measure of covariance complexity
% and its application to the evaluation of multivariate linear models. \emph{Comm.
% Stat. Theory and Methods} 19: 221-278
% Bozdogan, H. (2000). Akaike's Information Criteria and recent developments in
% Information Complexity. \emph{J. Math. Psych.} 44: 62-91

\references{
Mallows, C. L. (1973) Some Comments on CP. \emph{Technometrics} 15: 661–675.

Bozdogan, H. and Haughton, D.M.A. (1998) Information complexity criteria for
regression models. \emph{Comp. Stat. & Data Analysis} 28: 51-76.

}

\seealso{
\code{\link{AIC}} and \code{\link{BIC}} in \pkg{stats}, \code{\link{AICc}}
}

% \examples{}


\keyword{models}
