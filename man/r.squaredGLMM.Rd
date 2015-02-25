\name{r.squaredGLMM}
\alias{r.squaredGLMM}

\encoding{utf-8}

\newcommand{\Rsqx}{\ifelse{latex}{\eqn{R_{#1}^{2}}{R_#1^2}}{\ifelse{html}{\enc{R_#1²}{R_#1^2}}{R_#1^2}}}
\newcommand{\Rsq}{\ifelse{latex}{\eqn{R^{2}}{R^2}}{\ifelse{html}{\enc{R²}{R^2}}{R^2}}}


\title{Pseudo-R-squared for Generalized Mixed-Effect models}
\description{
Calculate conditional and marginal coefficient of determination for 
Generalized mixed-effect models (\Rsqx{GLMM}).

}

\usage{
r.squaredGLMM(x)
}

\arguments{
  \item{x}{a fitted linear model object. }
}

\value{
	\code{r.squaredGLMM} returns a numeric vector with two values for marginal
	and conditional \Rsqx{GLMM}.
}

\details{
For mixed-effects models, \Rsq can be categorized into two types. 
\bold{Marginal \Rsqx{GLMM}} represents the variance explained by fixed 
factors, and is defined as:

\deqn{
R_{GLMM(m)}^{2}= \frac{\sigma_f^2}{\sigma_f^2 + \sum_{l=1}^{u}\sigma_{l}^{2} + 
\sigma_e^2 + \sigma_d^2}
}{
R_GLMM(m)^2 = (\sigma_f^2) / (\sigma_f^2 + sum(\sigma_l^2) + \sigma_e^2 + 
\sigma_d^2
}

\bold{Conditional \Rsqx{GLMM}} is interpreted as variance explained by both
fixed and random factors (i.e. the entire model), and is calculated according
to the equation:

\deqn{
R_{GLMM(c)}^{2}= \frac{\sigma_f^2 + \sum_{l=1}^{u}\sigma_{l}^{2}}{\sigma_f^2 +
\sum_{l=1}^{u}\sigma_{l}^{2} + \sigma_e^2 + \sigma_d^2}
}{
R_GLMM(c)^2= (\sigma_f^2 + sum(\sigma_l^2)) / (\sigma_f^2 + sum(\sigma_l^2) +
\sigma_e^2 + \sigma_d^2 }

where \eqn{\sigma_f^2} is the variance of the fixed effect components, and
\eqn{\sum \sigma_{l}^{2}}{sum(\sigma_l^2)} is the sum of all
\ifelse{latex}{\eqn{u}}{} variance components (group, individual, etc.), 
\eqn{\sigma_l^2} is the variance due to additive dispersion and \eqn{\sigma_d^2}
is the distribution-specific variance. 

}


\note{
\Rsqx{GLMM} can be calculated also for fixed-effect models. In
the simpliest case of \acronym{OLS} it reduces to \code{var(fitted) /
(var(fitted) + deviance / 2)}. Unlike likelihood-ratio based \Rsq for 
\acronym{OLS}, value of this statistic differs from that of
the classical \Rsq.

Currently methods exist for classes: \code{mer}(\code{Mod}), \code{lme},
\code{glmmML} and (\code{g})\code{lm}.

See note in \code{\link{r.squaredLR}} help page for comment on using \Rsq in 
model selection. 

}

\author{
This implementation is based on \R code from \sQuote{Supporting Information} for
Nakagawa & Schielzeth (2012), and its extension by Paul Johnson.
}

\references{
Nakagawa, S, Schielzeth, H. (2013). A general and simple method for obtaining
\Rsq from Generalized Linear Mixed-effects Models. \emph{Methods in
Ecology and Evolution} 4: 133–142

Johnson, P.C.D. (2014) Extension Nakagawa & Schielzeth’s \Rsqx{GLMM} to random 
slopes models. \emph{Methods in Ecology and Evolution} 5: 44-946.
}

\seealso{
\code{\link{summary.lm}}, \code{\link{r.squaredLR}}
}

\examples{

\dontshow{ if(require(nlme)) \{   }
data(Orthodont, package = "nlme")

fm1 <- lme(distance ~ Sex * age, ~ 1 | Subject, data = Orthodont)

r.squaredGLMM(fm1)
r.squaredLR(fm1)
r.squaredLR(fm1, null.RE = TRUE)
\dontshow{  \}  }
}

\keyword{models}
