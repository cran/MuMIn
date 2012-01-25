\name{Model utilities}
\alias{MuMIn-model-utils}
\alias{coeffs}
\alias{model.names}
\alias{beta.weights}
\alias{tTable}
\alias{coefTable}
\alias{coefTable.default}
\alias{coefTable.lme}
\alias{getAllTerms}
\alias{getAllTerms.terms}

\encoding{utf-8}
\newcommand{\href}{#2}

\title{Model utility functions}
\description{
These functions extract or calculate various values from provided fitted model
	objects(s). They are mainly meant for internal use, but may be also useful
	for end-users.

\code{beta.weights} computes standardized coefficients (beta weights) for
	a model;

\code{coeffs} extracts model coefficients;

\code{getAllTerms} extracts independent variable names from a model object;

\code{coefTable} extracts a table of coefficients, standard errors and 
	associated degrees of freedom when possible;

\code{model.names} generates shorthand numeric names for one or several fitted
models.

}

\usage{
beta.weights(model)

coeffs(model)

getAllTerms(x, ...)
\method{getAllTerms}{terms}(x, offset = TRUE, intercept = FALSE, ...)

coefTable(model, ...)
\method{coefTable}{lme}(model, adjustSigma, ...)

model.names(object, ..., labels = NULL)
}
\arguments{
  \item{model}{a fitted model object. }

  \item{object}{a fitted model object or a list of such objects. }
  \item{x}{a fitted model object or a \code{formula}. }
  \item{offset}{should \sQuote{offset} terms be included? }
  \item{intercept}{should terms names include the intercept?}
  \item{labels}{optionally, a character vector with names of all the terms,
	e.g. from a global model. \code{model.names} enumerates the model terms in
	order of their appearance in the list and in the models. So, changing the
	order of the models would lead to different names. The argument 'labels' can
	be used to prevent this happening.
 }
  \item{\dots}{ For \code{model.names}, more fitted model objects.
	For \code{coefTable} arguments that are passed to appropriate \code{summary}
	method (e.g. \code{dispersion} parameter for \code{glm} may be used here).
	In other functions often not used.
  }
  \item{adjustSigma}{ See \code{\link[nlme]{summary.lme}}. }
}

\details{
	The functions \code{coeffs}, \code{getAllTerms} and \code{coefTable} provide
	interface between the model object and \code{model.avg} (and
	\code{dredge}). Custom methods can be written to provide support for
	additional classes of models. The vignette
	\sQuote{\ifelse{html}{\href{/library/MuMIn/doc/gamm.pdf}{Extending
	\pkg{MuMIn}'s functionality}}{Extending \pkg{MuMIn}'s functionality}}
	describes it in more detail.
}

\note{

\code{coeffs}'s value is in most cases identical to that returned by
\code{\link{coef}}, the only difference being it returns fixed
effects' coefficients for mixed models.

Use of \code{tTable} is deprecated in favour of \code{coefTable}.

}

\author{Kamil Barto\enc{ń}{n}}


\keyword{models}
\keyword{manip}
