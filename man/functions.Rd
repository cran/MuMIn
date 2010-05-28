\name{dRedging.functions}
\alias{coeffs}
\alias{beta.weights}
\alias{tTable}
\alias{getAllTerms}

\alias{coeffs.default}
\alias{coeffs.glmer}
\alias{coeffs.gls}
\alias{coeffs.lme}
\alias{coeffs.lmer}
\alias{coeffs.mer}
\alias{coeffs.spautolm}

\alias{getAllTerms.default}
\alias{getAllTerms.formula}
\alias{getAllTerms.glmer}
\alias{getAllTerms.lme}
\alias{getAllTerms.lmer}
\alias{getAllTerms.mer}

\alias{tTable.default}
\alias{tTable.gam}
\alias{tTable.glmer}
\alias{tTable.gls}
\alias{tTable.lme}
\alias{tTable.lmer}
\alias{tTable.mer}
\alias{tTable.spautolm}

\encoding{utf-8}

\title{Helper functions for package dRedging}
\description{
\code{beta.weights} - computes standardized coefficients (beta weights) for the model; \code{coeffs} - gets model coefficients; \code{getAllTerms} - extracts independent variables' names from a fitted model;
\code{tTable} - returns \code{\link[base]{summary}} table of a fitted model with coefficients, standard errors, and p-values.
}

\usage{
beta.weights(model) 
coeffs(model)
getAllTerms(x, ...)
tTable(model, ...)
}

\arguments{
  \item{model}{a fitted model object}
  \item{x}{a fitted model object or a \code{\link{formula}}}
  \item{\dots}{other arguments, not used}
}


\author{Kamil Bartoń}


\seealso{
\code{\link{dredge}}
}

\keyword{misc}
