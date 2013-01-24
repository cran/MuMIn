\name{r.squaredGLMM}
\alias{r.squaredGLMM}

\title{Pseudo-R-squared for Generalized Mixed-Effect models}
\description{
Calculate a conditional and marginal coefficient of determination for 
Generalized mixed-effect models (R\eqn{_{GLMM}^{2}}{_GLMM^2}).

}

\usage{
r.squaredGLMM(x, nullfx = NULL)
}

\arguments{
  \item{x}{a fitted linear model object. }
  \item{nullfx}{optionally, a fitted \emph{null} model including only intercept and 
	all the random effects of the reference model. }
}

\value{
	\code{r.squaredGLMM} returns a numeric vector with two values for marginal and 
		conditional R\eqn{_{GLMM}^{2}}{_GLMM^2}.
}

\details{
For mixed-effects models, R\eqn{^{2}}{^2} can be categorized into two types:
marginal and conditional. \bold{Marginal R\eqn{^{2}}{^2}} represents the variance explained by fixed 
factors, and is defined as:

\deqn{
R_{GLMM(m)}^{2}= \frac{\sigma_f^2}{\sigma_f^2
+ \sum_{l=1}^{u}\sigma_{l}^{2} + \sigma_\epsilon^2}
}{
R_GLMM(m)^2 = (\sigma_f^2) / (\sigma_f^2 + sum(\sigma_l^2) + \sigma_\epsilon^2
}

\bold{Conditional R\eqn{^{2}}{^2}} is interpreted as variance explained by both fixed and 
random factors (i.e . the entire model), and is calculated according to the equation:

\deqn{
R_{GLMM(c)}^{2}= \frac{\sigma_f^2 + \sum_{l=1}^{u}\sigma_{l}^{2}}{\sigma_f^2
+ \sum_{l=1}^{u}\sigma_{l}^{2} + \sigma_\epsilon^2}
}{
R_GLMM(c)^2= (\sigma_f^2 + sum(\sigma_l^2)) / (\sigma_f^2 + sum(\sigma_l^2) + \sigma_\epsilon^2
}

where \eqn{\sigma_f^2} is the variance of the fixed effect components, 
and \eqn{\sum \sigma_{l}^{2}}{sum(\sigma_l^2)} is the sum of all \ifelse{latex}{\eqn{u}}{} variance 
components (group, individual, etc.), and \eqn{\sigma_\epsilon^2} is the residual variance. 

}

\note{
R\eqn{_{GLMM}^{2}}{_GLMM^2} can be calculated also for fixed-effect models. In the 
simpliest case of \acronym{OLS} it reduces to \code{var(fitted) / (var(fitted) + 
deviance / 2)}.
Yet, unlike likelihood-ratio based R\eqn{^{2}}{^2} for \acronym{OLS}, value of this
statistic differs from that of the classical R\eqn{^{2}}{^2}.

Currently methods exist for classes: \code{mer}(\code{Mod}), \code{lme},
\code{glmmML} and (\code{g})\code{lm}.

See note in \code{\link{r.squaredLR}} help page for comment on using R\eqn{^{2}}{^2} 
in model selection. 

This implementation is based on \R code from \sQuote{Supporting Information} for Nakagawa & 
Schielzeth (2012).

}


\references{

Nakagawa, S, Schielzeth, H. (2012). A general and simple method for obtaining R\eqn{^{2}}{^2} from 
Generalized Linear Mixed-effects Models. \emph{Methods in Ecology and Evolution}: (online)
doi:10.1111/j.2041-210x.2012.00261.x

}

\seealso{
\code{\link{summary.lm}}, \code{\link{r.squaredLR}}
}

\examples{
library(lme4)
data(Orthodont, package = "nlme")

fm1 <- lmer(distance ~ Sex * age + (1 | Subject), data = Orthodont)

r.squaredGLMM(fm1)
r.squaredLR(fm1)
r.squaredLR(fm1, null.RE = TRUE)

}

\keyword{models}
