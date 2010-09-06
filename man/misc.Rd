\name{miscellaneous}
\alias{Miscellaneous}

\alias{coeffs}
\alias{beta.weights}
\alias{tTable}
\alias{getAllTerms}

%%\alias{coeffs.default}
%%\alias{coeffs.glmer}
%%\alias{coeffs.gls}
%%\alias{coeffs.lme}
%%\alias{coeffs.lmer}
%%\alias{coeffs.mer}
%%\alias{coeffs.spautolm}

%%\alias{getAllTerms.default}
%%\alias{getAllTerms.formula}
\alias{getAllTerms.terms}
%%\alias{getAllTerms.glmer}
%%\alias{getAllTerms.lme}
%%\alias{getAllTerms.lmer}
%%\alias{getAllTerms.mer}

%%\alias{tTable.default}
%%\alias{tTable.gam}
%%\alias{tTable.glmer}
%%\alias{tTable.gls}
%%\alias{tTable.lme}
%%\alias{tTable.lmer}
%%\alias{tTable.mer}
%%\alias{tTable.spautolm}

\alias{Weights}
\alias{cbindDataFrameList}
\alias{rbindDataFrameList}

\encoding{utf-8}

\title{Helper functions}
\description{

\code{beta.weights} - computes standardized coefficients (beta weights) for
	a model;

\code{coeffs} - extracts model coefficients;

\code{getAllTerms} - extracts independent variable names from a model object;

\code{tTable} - extracts a table of coefficients,
	standard errors, and p-values from a model object;

\code{Weights} - calculates Akaike weights (normalized relative likelihoods)

}

\usage{
beta.weights(model)
coeffs(model)
getAllTerms(x, ...)
\method{getAllTerms}{terms}(x, offset = TRUE, ...)
tTable(model, ...)
Weights(aic, ...)

cbindDataFrameList(x)
rbindDataFrameList(x)

}

\arguments{
  \item{model}{a fitted model object}
  \item{x}{a fitted model object or a \code{\link{formula}}.
  for \code{*bindDataFrameList}, a list of \code{data.frames}
  }
  \item{offset}{should \sQuote{offset} terms be included?}
  \item{\dots}{other arguments, not used}
  \item{aic}{a vector of AIC (or other information criterion) values }
}

\details{
	The functions \code{coeffs}, \code{getAllTerms} and \code{tTable} provide an
	interface between the model and \code{model.avg} (as well as \code{dredge}).
	Custom methods can be written to provide support for additional classes of
	models. Also, a \code{logLik} method must exist for object.
}

\note{

\code{coeffs}'s value is in most cases identical to that returned by
\code{\link{coef}}, the only difference is that it returns fixed
effects' coefficients for mixed models.

Functions \code{*bindDataFrameList} are not exported from the name space,
use \code{MuMIn:::cbindDataFrameList} to access them.

}


\author{Kamil Bartoń}


\seealso{
\code{\link{dredge}}
}


\keyword{models}
\keyword{manip}
