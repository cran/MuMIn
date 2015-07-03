\name{subset.model.selection}
\alias{subset.model.selection}
\alias{[.model.selection}
\alias{[[.model.selection}
\alias{has}

\newcommand{\bq}{\verb{`}\code{#1}\verb{`}}

\encoding{utf-8}
\title{Subsetting model selection table}
\description{Extract subset of a model selection table.
}

\usage{
\method{subset}{model.selection}(x, subset, select, recalc.weights = TRUE, recalc.delta = FALSE, ...)
\method{[}{model.selection}(x, i, j, recalc.weights = TRUE, recalc.delta = FALSE, ...)
\method{[[}{model.selection}(x, ..., exact = TRUE) 

}

\arguments{
  \item{x}{a \code{model.selection} object to be subsetted.}
  \item{subset,select}{logical expressions indicating columns and rows to keep.
    See \code{\link{subset}}. }
  \item{i,j}{indices specifying elements to extract. }
  \item{recalc.weights}{logical value specyfying whether Akaike weights
        should be normalized across the new set of models to sum to one.
    }
    \item{recalc.delta}{logical value specyfying whether 
        \ifelse{latex}{\eqn{\Delta_{IC}}}{\enc{Δ}{Delta}_IC}
        should be calculated for the new set of models (not done by default).
    }
    \item{exact}{logical, see \code{\link[=Extract]{[}}}.
  \item{\dots}{further arguments passed to \code{\link{[.data.frame}} 
    (\code{drop}). }
}

\value{
A \code{model.selection} object containing only the selected models (rows).
If columns are selected (via argument \code{select} or second index 
\code{x[, j]}) and not all if not all essential columns (i.e. all except 
"varying" and "extra") are present in the result, a plain \code{data.frame} is 
returned. Similarly, modifying values in the essential columns with \code{[<-}, 
\code{[[<-} or \code{$<-} produces a regular data frame. 

}

\details{
Unlike the method for \code{data.frame}, single bracket extraction with only 
one index \code{x[i]} selects rows (models) rather than columns.

To select rows according to presence or absence of the variables (rather than
their value), a pseudo-function \code{has} may be used with \code{subset}, e.g.
\code{subset(x, has(a, !b))}
will select rows  with \emph{a} \bold{and} without \emph{b} (this is
equivalent to \code{!is.na(a) & is.na(b)}). \code{has} can take any number of
arguments.

Complex model terms need to be enclosed within curly brackets 
 (e.g \code{{s(a,k=2)}}), except for within \code{has}. Backticks-quoting is 
 also possible, but then the name must match exactly (including whitespace) 
 the term name as returned by \code{getAllTerms}.

To select rows where one variable can be present conditional on the presence of 
other variable(s), the function \code{dc} (\bold{d}ependency \bold{c}hain) can 
be used. 
\code{dc} takes any number of variables as arguments, and allows a variable to be 
included only if all the preceding arguments are also included (e.g. \code{subset = 
dc(a, b, c)} allows for models of form \code{a}, \code{a+b} and \code{a+b+c} but not 
\code{b}, \code{c}, \code{b+c} or \code{a+c}).

}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{dredge}}, \code{\link{subset}} and \code{\link{[.data.frame}} for
subsetting and extracting from \code{data.frame}s.
}

\examples{
fm1 <- lm(formula = y ~ X1 + X2 + X3 + X4, data = Cement, na.action = na.fail)

# generate models where each variable is included only if the previous
# are included too, e.g. X2 only if X1 is there, and X3 only if X2 and X1
dredge(fm1, subset = dc(X1, X2, X3, X4))

# which is equivalent to
# dredge(fm1, subset = (!X2 | X1) & (!X3 | X2) & (!X4 | X3))

# alternatively, generate "all possible" combinations
ms0 <- dredge(fm1)
# ...and afterwards select the subset of models
subset(ms0, dc(X1, X2, X3, X4))
# which is equivalent to
# subset(ms0, (has(!X2) | has(X1)) & (has(!X3) | has(X2)) & (has(!X4) | has(X3)))

# Different ways of finding a confidence set of models:
# delta(AIC) cutoff
subset(ms0, delta <= 4, recalc.weights = FALSE)
# cumulative sum of Akaike weights
subset(ms0, cumsum(weight) <= .95, recalc.weights = FALSE)
# relative likelihood
subset(ms0, (weight / weight[1]) > (1/8), recalc.weights = FALSE)
}


\keyword{manip}
