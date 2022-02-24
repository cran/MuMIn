\encoding{utf-8}
\name{plot.model.selection}
\alias{plot.model.selection}
\title{Visualize model selection table}
\usage{
\method{plot}{model.selection}(x, ylab = NULL, xlab = NULL, labels = attr(x, "terms"),
  labAsExpr = TRUE, mar.adj = TRUE,
  col = c("SlateGray", "SlateGray2"), col2 = "white", border = par("col"),
  par.lab = NULL, par.vlab = NULL,
  axes = TRUE, ann = TRUE, ...)
}
\arguments{
\item{x}{a \code{"model.selection"} object.}

\item{xlab, ylab}{labels for the \emph{x} and \emph{y} axes.}

\item{labels}{optional, a character vector or an expression containing model
term labels (to appear on top side of the plot).
Its length must be equal to number of model terms in the table. Defaults
to the model term names.}

\item{labAsExpr}{a logical indicating whether the character labels should be
interpreted (\code{\link{parse}}d) as \R expressions.}

\item{mar.adj}{logical indicating whether the upper margin should be
expanded to fit the labels}

\item{col, col2}{vector of colors for columns (if more than one \code{col} is
given, columns will be filled with alternating colors). If \code{col2} is
specified cells will be filled with gradient from \code{col} to \code{col2}.
\code{col2 = NA} produces no gradient.}

\item{border}{border color for cells and axes.}

\item{par.lab, par.vlab}{optional lists of \link[=par]{graphical parameters} for the
term labels (top axis) and model names (right axis), respectively.}

\item{axes, ann}{logical values indicating whether the axis and annotation
should appear on the plot.}

\item{\dots}{further \link[=par]{graphical parameters} to be set for the plot.}
}
\value{
A numeric vector of model weights.
}
\description{
Produces a graphical representation of model weights per model term.
}
\examples{
# For examples, see \sQuote{[MuMIn-package]}
}
\seealso{
\code{\link{plot.default}}, \code{\link{par}}
}
\author{
Kamil Barto\enc{ń}{n}
}
\keyword{hplot}
