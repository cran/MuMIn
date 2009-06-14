\name{AICc}
\alias{AICc}
\encoding{utf-8}
\title{Second-order Akaike Information Criterion}
\description{
Calculates second-order Akaike information criterion for one or several fitted model objects (AIC for small samples).
}

\usage{
AICc(object, ..., k = 2)
}

\arguments{
  \item{object}{a fitted model object}
  \item{\dots}{optionally more fitted model objects}
  \item{k}{the \dQuote{penalty} per parameter to be used; the default \code{k} = 2 is the classical \code{\link[stats]{AIC}}}
}


\value{
	If just one object is provided, returns a numeric value with the corresponding AICc; 
	if more than one object are provided, returns a data.frame with rows corresponding to the objects and columns representing the number of parameters in the model (df), AICc and the \code{\link[stats]{AIC}}.

}
\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel inference: a practical information-theoretic approach}. 2nd ed. 
}

\author{ Kamil Bartoń}

\seealso{ 
Akaike's An Information Criterion: \code{\link[stats]{AIC}}
}

\keyword{models}