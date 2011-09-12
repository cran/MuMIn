\name{par.avg}
\alias{par.avg}
\encoding{utf-8}

\title{Parameter averaging}
\description{
Averages single model coefficient based on provided weights
}

\usage{
par.avg(x, se, weight, df = NULL, alpha = 0.05, 
    revised.var = TRUE)
}

\arguments{
  \item{x}{vector of parameters}
  \item{se}{vector of standard errors}
  \item{weight}{vector of weights}
  \item{df}{(optional) vector of degrees of freedom}
  \item{alpha}{significance level for calculating confidence intervals}
  \item{revised.var}{logical, should the revised formula for standard errors be
	used? See \sQuote{Details}.}
}

\value{
  \code{par.avg} returns a vector with named elements:
  \item{Coefficient}{model coefficients}
  \item{SE}{unconditional standard error}
  \item{Adjusted SE}{adjusted standard error}
  \item{Lower CI, Upper CI}{unconditional confidence intervals}
}

\details{
Unconditional standard errors are square root of the variance estimator,
calculated either according to the original formula in
Burnham and Anderson (2002, p. 160, equation 4.7),
or a newer, revised formula from Burnham and Anderson (2004, equation 4)
(if \code{revised.var = TRUE}, this is the default).
If degrees of freedom are given, the confidence intervals are based on adjusted
standard error estimator (Burnham and Anderson 2002, page 164).
}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed.

Burnham, K. P. and Anderson, D. R. (2004). \emph{Multimodel inference -
understanding AIC and BIC in model selection}. Sociological Methods & Research
33(2): 261-304.
}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{model.avg}} for model averaging.
}


\keyword{models}
