\name{merge.model.selection}
\alias{merge.model.selection}
\alias{rbind.model.selection}
\alias{append.model.selection}

\encoding{utf-8}
\title{Combine model selection tables}
\description{Combine two or more model selection tables. }

\usage{
\method{merge}{model.selection}(x, y, suffixes = c(".x", ".y"), ...)
% \method{rbind}{model.selection}(..., deparse.level = 1)
}

\arguments{
  \item{x, y}{\code{model.selection} objects to be combined.}
  \item{suffixes}{a character vector with two elements that are appended
	respectively to row names of the combined tables}
  \item{\dots}{ignored}
% \item{\dots}{\code{model.selection} objects for \code{rbind}, ignored in 
% \code{merge}}
%  \item{deparse.level}{currently ignored.}
}

\value{
A \code{model.selection} object containing the selected models (rows).
}

\note{
Both \ifelse{latex}{\eqn{\Delta_{IC}}}{\enc{Δ}{Delta}_IC} values and 
\emph{Akaike weights} are recalculated in the resulting tables.

Models in the combined model selection tables must be comparable, i.e. fitted to
	the same data, however only very basic checking is done to verify that. The 
	models must also be ranked by the same information criterion.
	
Unlike the \code{merge} method for \code{data.frame}, this method appends
	second table to the first (similarly to \code{rbind}).
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{dredge}}, \code{\link{model.sel}}, \code{\link{merge}}, 
	\code{\link{rbind}}.
}

\examples{
\dontrun{

require(mgcv)
data(Beetle)

ms1 <- dredge(glm(Prop ~ dose + I(dose^2) + log(dose) + I(log(dose)^2),
    data = Beetle, family = binomial, na.action = na.fail))
	
fm2 <- gam(Prop ~ s(dose, k = 3), data = Beetle, family = binomial)

merge(ms1, model.sel(fm2))
}
}


\keyword{manip}
