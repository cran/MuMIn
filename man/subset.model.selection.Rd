\name{subset.model.selection}
\alias{subset.model.selection}
\alias{[.model.selection}
\alias{has}

\encoding{utf-8}
\title{Subsetting model selection table}
\description{Return subsets of a model selection table returned by
	\code{dredge} or \code{model.sel}.
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

\note{
Unlike the method for \code{data.frame}, extracting with only one index
(i.e. \code{x[i]}) will select rows rather than columns.

To select rows according to presence or absence of the variables (rather than
their value), a pseudo-function \code{has} may be used, e.g.
\code{subset(x, has(a, !b))}
will select rows  with \emph{a} \bold{and} without \emph{b} (this is
equivalent to \code{!is.na(a) & is.na(b)}). \code{has} can take any number of
arguments. Importantly, the \code{has()} notation cannot be used in the
\code{subset} argument for \code{dredge}, where the variable names should be
given directly, with the same effect.
%\code{subset = has(a, !b)}
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{dredge}}, \code{\link{subset}} and \code{\link{[.data.frame}} for
subsetting and extracting from \code{data.frame}s.
}

\examples{
data(Cement)
fm1 <- lm(formula = y ~ X1 + X2 + X3 + X4, data = Cement)

# generate models where each variable is included only if the previous
# are included too, e.g. X2 only if X1 is there, and X3 only if X2 and X1
dredge(fm1, subset = (!X2 | X1) & (!X3 | X2) & (!X4 | X3))

# alternatively, generate "all possible" combinations
ms0 <- dredge(fm1)
# ...and afterwards select the subset of models
subset(ms0, (has(!X2) | has(X1)) & (has(!X3) | has(X2)) & (has(!X4) | has(X3)))
\dontrun{
# this way the expression may be more clear
subset(ms0, has(X1, X2, X3, X4) | has(X1, X2, X3) | has(X1, X2) | has(X1)
    | (df == 2))
}

# Different ways of finding a confidence set of models:
# delta(AIC) cutoff
subset(ms0, delta <= 4, recalc.weights = FALSE)
# cumulative sum of Akaike weights
subset(ms0, cumsum(weight) <= .95, recalc.weights = FALSE)
# relative likelihood
subset(ms0, (weight / weight[1]) > (1/8), recalc.weights = FALSE)
}


\keyword{manip}
