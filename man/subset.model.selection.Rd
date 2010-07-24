\name{subset.model.selection}
\alias{subset.model.selection}
\alias{[.model.selection}

\encoding{utf-8}
\title{Subsetting model selection table}
\description{Return subsets of a model selection table returned by
\code{dredge}.
}

\usage{

\method{subset}{model.selection}(x, subset, select, recalc.weights = TRUE, ...)
\method{[}{model.selection}(x, i, j, recalc.weights = TRUE, ...)

}

\arguments{
  \item{x}{a \code{model.selection} object to be subsetted.}
  \item{subset,select}{logical expressions indicating columns and rows to keep.
	See \code{\link{subset}}. }
  \item{i,j}{indices specifying elements to extract. }
  \item{recalc.weights}{logical value specyfying whether Akaike weights should
	be normalized across the new set of models to sum to one.
	}
  \item{\dots}{further arguments passed to  \code{\link{[.data.frame}}. }
}


\value{
A \code{model.selection} object containing only the selected models (rows).
When columns are selected (arguments \code{select} or \code{j} are provided), a
plain \code{data.frame} is returned.
}

\author{Kamil Bartoń}


\seealso{
\code{\link{dredge}}, \code{\link{subset}} and \code{\link{[.data.frame}} for
subsetting and extracting from data.frames.
}


\keyword{manip}
