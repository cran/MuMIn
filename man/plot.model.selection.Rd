\encoding{utf-8}
\name{plot.model.selection}
\alias{plot.model.selection}
\title{Visualize model selection table}
\usage{
\method{plot}{model.selection}(
  x,
  ylab = NULL, xlab = NULL, main = "Model selection table",
  labels = NULL, terms = NULL, labAsExpr = TRUE,
  vlabels = rownames(x), mar.adj = TRUE,
  col = NULL, col.mode = 2,
  bg = "white", border = par("col"),
  par.lab = NULL, par.vlab = NULL,
  axes = TRUE, ann = TRUE,
  ...
)
}
\arguments{
\item{x}{a \code{"model.selection"} object.}

\item{xlab, ylab, main}{labels for the \emph{x} and \emph{y} axes, and the
main title for the plot.}

\item{labels}{optional, a character vector or an expression containing model
term labels (to appear on top side of the plot).
Its length must be equal to number of displayed model terms. Defaults
to the model term names.}

\item{terms}{which terms to include (default \code{NULL} means all terms).}

\item{labAsExpr}{logical, indicating whether the term names should be
interpreted (\lcode{parse}d) as \R expressions for prettier labels.
See also \lcode{plotmath}.}

\item{vlabels}{alternative labels for the table rows (i.e. model names)}

\item{mar.adj}{logical indicating whether the top and right margin should be
enlarged if necessary to fit the labels.}

\item{col}{vector or a \code{matrix} of colours for the non-empty grid cells.
See 'Details'.
If \code{col} is given as a matrix, the colours are applied to rows and columns. How it
is done is governed by the argument \code{col.mode}.}

\item{col.mode}{either numeric or \code{"value"}, specifies cell colouring mode.
See 'Details'.}

\item{bg}{background colour for the empty cells.}

\item{border}{border colour for cells and axes.}

\item{par.lab, par.vlab}{optional lists of arguments and
\link[=par]{graphical parameters} for drawing
term labels (top axis) and model names (right axis), respectively.
Items of \code{par.lab} are passed as arguments to \lcode{mtext}, and those of
\code{par.vlab} are passed to \lcode{axis}.}

\item{axes, ann}{logical values indicating whether the axis and annotation
should appear on the plot.}

\item{\dots}{further \link[=par]{graphical parameters} to be set for the plot.}
}
\description{
Produces a graphical representation of model weights and terms.
}
\details{
\subsection{Colours}{

If \code{col.mode = 0}, the colours are recycled: if \code{col} is a matrix,
recycling takes place both per row and per column. If \code{col.mode > 0}, the
colour values in the columns are interpolated and assigned according to
the model weights. Higher values shift the colours for models with lower
model weights more forward. See also \lcode{colorRamp}. If \code{col.mode < 0} or
\code{"value"} (partially matched, case-insensitive) and \code{col} has two or more
elements, colours are used to represent coefficient values: the first
element in \code{col} is used for categorical predictors, the rest for
continuous values. The default is \code{grey} for factors and
\link[=hcl.colors]{HCL palette} \code{"Blue-Red 3"} otherwise, ranging from blue
for negative values to red for positive ones.
}

The following arguments are useful for adjusting label size and
position in \code{par.lab} and \code{par.vlab} : \code{cex}, \code{las} (see \lcode{par}),
\code{line} and \code{hadj} (see \lcode{mtext} and \lcode{axis}).
}
\examples{
\dontshow{if (require(graphics)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
ms <- dredge(lm(formula = y ~ ., data = Cement, na.action = na.fail))
plot(ms,
    # colours by coefficient value:
    col.mode = "value",
    par.lab = list(las = 2, line = 1.2, cex = 1),
    bg = "gray30",
    # change labels for the models to Akaike weights:
    vlabels = parse(text = paste("omega ==", round(Weights(ms), 2)))
    )
plot(ms, col = 2:3, col.mode = 0) # colour recycled by row
plot(ms, col = cbind(2:3, 4:5), col.mode = 0) # colour recycled by row and column
plot(ms, col = 2:3, col.mode = 1) # colour gradient by model weight
\dontshow{\}) # examplesIf}
}
\seealso{
\lcode{plot.default}, \lcode{par}, \code{\link{MuMIn-package}}
}
\author{
Kamil Barto\enc{ń}{n}
}
\keyword{hplot}
