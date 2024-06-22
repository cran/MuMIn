\name{r.squaredLR}
\alias{r.squaredLR}
\alias{null.fit}
\encoding{utf-8}

%% \newcommand{\mydequation}{\ifelse{latex}{\deqn{#1}{#3}}{\ifelse{html}{\enc{\deqn{#1}}{#3}}{#3}}}
%% \newcommand{\myequation}{\ifelse{latex}{\eqn{#1}{#3}}{\ifelse{html}{\enc{\eqn{#1}}{#3}}{#3}}}

%% \newcommand{\logLik}{\myequation{\log\mathcal{L}(#1)}{XXXX}{logL(#1)}}
%% \newcommand{\Rsq}{\myequation{R^{2}}{XXXX}{R^2}}
%% \newcommand{\Rsqx}{\myequation{R_{#1}^{2}}{XXXX}{R_#1^2}}



\title{Likelihood-ratio based pseudo-R-squared}
\description{
Calculate a coefficient of determination based on the likelihood-ratio test
(\Rsqx{LR}).
}

\usage{
r.squaredLR(object, null = NULL, null.RE = FALSE, ...)

null.fit(object, evaluate = FALSE, RE.keep = FALSE, envir = NULL, ...)
}

\arguments{
  \item{object}{a fitted model object. }
  \item{null}{a fitted \emph{null} model. If not provided, \code{null.fit} will
	be used to construct it. \code{null.fit}'s capabilities are limited to only 
	a few model classes, for others the \emph{null} model has to be specified
	manually. }
  \item{null.RE}{ logical, should the null model contain random factors?  Only 
    used if no \emph{null} model is given, otherwise omitted, with a warning. }
  \item{evaluate}{if \code{TRUE} evaluate the fitted model object else return
	the call. }
  \item{RE.keep}{if \code{TRUE}, the random effects of the original model are 
	included. }
  \item{envir}{the environment in which the \emph{null} model is to be
    evaluated, defaults to the environment of the original model's formula. }
  \item{\dots}{further arguments, of which only \code{x} would be used, to 
    maintain compatibility with older versions (\code{x} has been replaced 
	with \code{object}).}
}

\value{
	\code{r.squaredLR} returns a value of \Rsqx{LR}, and the
	attribute \code{"adj.r.squared"} gives the Nagelkerke's modified statistic.
	Note that this is not the same as nor equivalent to the classical 
	\sQuote{adjusted R squared}.

	\code{null.fit} returns the fitted \emph{null} model object (if
	\code{evaluate = TRUE}) or an unevaluated call to fit a \emph{null} model.

}

\details{

This statistic is is one of the several proposed pseudo-\Rsq's for
nonlinear regression models. It is based on an improvement from \emph{null}
(intercept only) model to the fitted model, and calculated as

\mydequation{
  R_{LR}^{2}=1-\exp(-\frac{2}{n}(\log\mathcal{L}(x)-\log\mathcal{L}(0)))
}{R² = 1 - exp(-2/n * ㏒ℓ(x) - ㏒ℓ(0))
}{R^2 = 1 - exp(-2/n * logL(x) - logL(0))
}

where \logLik{x} and \logLik{0} are the log-likelihoods of the
fitted and the \emph{null} model respectively.
\acronym{ML} estimates are used if models have been 
fitted by \acronym{RE}stricted \acronym{ML} (by calling \code{logLik} with argument 
\code{REML = FALSE}). Note that the \emph{null} model can include the random 
factors of the original model, in which case the statistic represents the 
\sQuote{variance explained} by fixed effects.

For \acronym{OLS} models the value is consistent with classical \Rsq. In some 
cases (e.g. in logistic regression), the maximum \Rsqx{LR} is less than one. 
The modification proposed by Nagelkerke (1991) adjusts the \Rsqx{LR} to achieve 
1 at its maximum:
\myequation{\bar{R}^{2} = R_{LR}^{2} / \max(R_{LR}^{2})
}{R̅² = R² / max(R²)
}{Radj^2 = R^2 / max(R^2)
}
where
\myequation{\max(R_{LR}^{2}) = 1 - \exp(\frac{2}{n}\log\mathcal{L}(\textrm{0}))
}{max(R²) = 1 - exp(2 / n * ㏒ℓ(0))
}{max(R^2) = 1 - exp(2 / n * logL(0))
}.

\code{null.fit} tries to guess the \emph{null} model call, given the provided
fitted model object. This would be usually a \code{glm}. The function will give
an error for an unrecognised class.

}

\note{
\Rsq is a useful goodness-of-fit measure as it has the interpretation
of the proportion of the variance \sQuote{explained}, but it performs poorly in
model selection, and is not suitable for use in the same way as the information
criteria.
}


\references{

Cox, D. R. and Snell, E. J. 1989 \emph{The analysis of binary data}, 2nd ed.
London, Chapman and Hall.

Magee, L. 1990 \Rsq measures based on Wald and likelihood ratio joint
significance tests. \emph{Amer. Stat.} \bold{44}, 250--253.

Nagelkerke, N. J. D. 1991 A note on a general definition of the coefficient of
determination. \emph{Biometrika} \bold{78}, 691--692.
}

\seealso{
\code{\link{summary.lm}}, \code{\link{r.squaredGLMM}}

\code{\link[performance]{r2}} from package \CRANpkg{performance} calculates
many different types of \Rsq.

}

\keyword{models}
