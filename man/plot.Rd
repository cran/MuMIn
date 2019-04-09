\name{plot.model.selection}
\alias{plot.model.selection}
\encoding{utf-8}

\title{Visualize model selection table}
\description{
Produces a graphical representation of model weights per model term.
}

\usage{
\method{plot}{model.selection}(x, ylab = NULL, xlab = NULL,
    labels = attr(x, "terms"), labAsExpr = FALSE,
    col = c("SlateGray", "SlateGray2"), col2 = "white", border = par("col"),
    par.lab = NULL, par.vlab = NULL,
    axes = TRUE, ann = TRUE, ...) 
}
\arguments{
    \item{x}{a \code{"model.selection"} object. }
    \item{xlab, ylab}{labels for the x and y axis. }
	\item{labels}{optional, a character vector or an expression containing model
		term labels (to appear on top side of the plot).
		Its length must be equal to number of model terms in the table. Defaults
		to model term names. }
	\item{labAsExpr}{a logical indicating whether the character labels should be
		interpreted (\code{\link{parse}}d) as \R expressions. }
	\item{col, col2}{vector of colors for columns (if more than one \code{col}
		is given, columns will be filled with alternating colors). If
		\code{col2} is specified cells will be filled with gradient from \code{col}
		to \code{col2}. Set \code{col2} to \code{NA} for no gradient. }
    \item{border}{border color for cells and axes. }
	\item{par.lab, par.vlab}{ optional lists or parameters for term labels (top
		axis) and model names (right axis), respectively. }
	\item{axes, ann}{ logical values indicating whether the axis and annotation
		should appear on the plot.}
	\item{\dots}{further \code{\link[=par]{graphical parameters}} to be set for
		the plot (see \code{\link{par}}).}
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{plot.default}, \link{par}}

For examples, see \sQuote{\link{MuMIn-package}}
}

\keyword{hplot}
