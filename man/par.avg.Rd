\name{par.avg}
\alias{par.avg}
\encoding{utf-8}

\title{Parameter averaging}
\description{
Averages single parameter based on provided weights
}

\usage{
par.avg(x, se, npar, weight, alpha = 0.05)
}

\arguments{
  \item{x}{vector of parameters}
  \item{se}{vector of standard errors}
  \item{npar}{vector giving numbers of estimated parameters}
  \item{weight}{vector of weights}
  \item{alpha}{significance level for calculatinq confidence intervals}
}

\value{
  par.avg returns a vector with named elements:
  \item{Coefficient}{model coefficients}
  \item{Variance}{unconditional variance of coefficients}
  \item{Unconditional SE}{ }
  \item{Lower CI, Upper CI}{relative variable importances}
}

\references{
Burnham, K. P. and Anderson, D. R (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed. }

\author{Kamil Bartoń}

\seealso{
\code{\link{model.avg}} for averaging models.
}


\keyword{models}
