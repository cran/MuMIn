\name{Miscellaneous}
\alias{Miscellaneous}
\alias{miscellaneous}
\alias{coeffs}
\alias{model.names}
\alias{beta.weights}
\alias{tTable}
\alias{getAllTerms}
\alias{getAllTerms.terms}


\alias{cbindDataFrameList}
\alias{rbindDataFrameList}

\encoding{utf-8}

\title{Helper functions}
\description{

\code{beta.weights} computes standardized coefficients (beta weights) for
	a model;

\code{coeffs} extracts model coefficients;

\code{getAllTerms} extracts independent variable names from a model object;

\code{tTable} extracts a table of coefficients, standard errors, and p-values
from a model object;

\code{model.names} generates shorthand numeric names for one or several fitted
models.

}

\usage{
beta.weights(model)
coeffs(model)
getAllTerms(x, ...)
\method{getAllTerms}{terms}(x, offset = TRUE, intercept = FALSE, ...)
tTable(model, ...)
model.names(object, ..., labels = NULL)

cbindDataFrameList(x)
rbindDataFrameList(x)
}

\arguments{
  \item{model}{a fitted model object. }
  \item{object}{a fitted model object or a list of such objects. }
  \item{x}{a fitted model object or a \code{\link{formula}}.
	for \code{*bindDataFrameList}, a list of \code{data.frames}. }
  \item{offset}{should \sQuote{offset} terms be included? }
  \item{intercept}{should terms names include the intercept?}
  \item{labels}{optionally, a character vector with names of all the terms,
	e.g. from a global model. \code{model.names} enumerates the model terms in
	order of their appearance in the list and in the models. So, changing the
	order of the models would lead to different names. The argument 'labels' can
	be used to prevent this happening.
 }
  \item{\dots}{ For \code{model.names}, more fitted model objects.
	In other functions often not used.
  }
}


\details{
	The functions \code{coeffs}, \code{getAllTerms} and \code{tTable} provide
	interface between the fitted model object and \code{model.avg} (and
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

Whimsically, the functions \code{*bindDataFrameList} are not exported from the
name space, \cr use \code{MuMIn:::cbindDataFrameList} to access them.

}

\author{Kamil Barto\enc{ń}{n}}


\keyword{models}
\keyword{manip}
