\name{r.squaredLR}
\alias{r.squaredLR}
\alias{null.fit}

\title{Likelihood-ratio based pseudo-R-squared}
\description{
Calculate a coefficient of determination based on the likelihood-ratio test
(\ifelse{latex}{\eqn{R_{LR}^{2}}}{R_LR-squared}).
}

\usage{
r.squaredLR(x, null = null.fit(x, TRUE, parent.frame()))

null.fit(x, evaluate = FALSE, envir = environment(as.formula(formula(x))))
}

\arguments{
  \item{x}{a fitted model object. }
  \item{null}{a fitted \emph{null} model. If not provided, \code{null.fit} will
	be used to construct it. Its capabilities are limited to only a few model 
	classes, for others the \emph{null} model has to be specified manually. }
  \item{evaluate}{If \code{TRUE} evaluate the fitted model object else return
	the call. }
  \item{envir}{the environment in which the \emph{null} model is to be
  evaluated, defaults to the environment of the original model's formula. }
}

\value{
	\code{r.squaredLR} returns a value of \eqn{R_{LR}^{2}}{R_LR^2}, and the
	attribute \code{"adj.r.squared"} gives the Nagelkerke's modified statistic.
	Note that this is not the same as the classical \sQuote{adjusted R squared}.

	\code{null.fit} returns the fitted \emph{null} model object (if
	\code{evaluate = TRUE}) or an unevaluated call to fit a \emph{null} model.

}

\details{

This statistic is is one of the several proposed pseudo-R-squared's for
nonlinear regression models. It is based on an improvement from \emph{null} 
(intercept only) model to the fitted model, and calculated as
\deqn{
R_{LR}^{2}=1-\exp(-\frac{2}{n}(\log\mathit{Lik}(x)-\log\mathit{Lik}(0)))
}{R^2 = 1 - exp(-2/n * logLik(x) - logLik(0))}
where log\emph{Lik(x)} and log\emph{Lik(0)} are the log-likelihoods of the fitted
and the \emph{null} model respectively.

For OLS models the value is consistent with classical \eqn{R^{2}}{R-squared}.
In some cases (e.g. in logistic regression), the maximum \eqn{R_{LR}^{2}}{R_LR^2}
is less than one. The modification proposed by Nagelkerke (1991) adjusts the
\eqn{R_{LR}^{2}}{R^2}
to achieve 1 at its maximum:
\eqn{
\bar{R}^{2} = R_{LR}^{2} / \max(R_{LR}^{2})
}{
Radj^2 = R^2 / max(R^2)
}
where
\eqn{
max(R_{LR}^{2}) = 1 - \exp(\frac{2}{n}\log\mathit{Lik}(\textrm{0}))
}{max(R^2) = 1 - exp(2 / n * logLik(0))}.

\code{null.fit} tries to guess the \emph{null} model call,
given the provided fitted model object. This would be usually a \code{glm}. The 
function will give an error for an unrecognized class.

}

\note{
R\eqn{^{2}}{^2} is a useful goodness-of-fit measure as it has the
interpretation of the proportion of the variance \sQuote{explained}, but it
performs poorly in model selection, and is not suitable for use in the same way
as the information criterions.
}


\references{
Cox, D. R. and Snell, E. J. (1989) \emph{The analysis of binary data}, 2nd ed.
London, Chapman and Hall

Magee, L. (1990) R\eqn{^{2}}{^2} measures based on Wald and likelihood ratio joint
significance tests. \emph{Amer. Stat.} 44: 250-253

Nagelkerke, N. J. D. (1991) A note on a general definition of the coefficient of
determination. \emph{Biometrika} 78: 691-692
}

\seealso{
\code{\link{summary.lm}}
}

\keyword{models}
