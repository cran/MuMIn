\name{r.squaredGLMM}
\alias{r.squaredGLMM}
\alias{r.squaredGLMM.merMod}

\encoding{utf-8}

%% \newcommand{\CRANpkg}{\href{https://cran.r-project.org/package=#1}{\pkg{#1}}} 
%% \newcommand{\mydequation}{\ifelse{latex}{\deqn{#1}{#3}}{\ifelse{html}{\enc{\deqn{#1}}{#3}}{#3}}}
%% \newcommand{\myequation}{\ifelse{latex}{\eqn{#1}{#3}}{\ifelse{html}{\enc{\eqn{#1}}{#3}}{#3}}}
%% 
%% \newcommand{\logLik}{\myequation{\log\mathcal{L}(#1)}{XXXX}{logL(#1)}}
%% \newcommand{\Rsq}{\myequation{R^{2}}{XXXX}{R^2}}
%% \newcommand{\Rsqx}{\myequation{R_{#1}^{2}}{XXXX}{R_#1^2}}

\title{Pseudo-R-squared for Generalized Mixed-Effect models}
\description{
Calculate conditional and marginal coefficient of determination for 
Generalized mixed-effect models (\Rsqx{GLMM}).
}

\usage{
r.squaredGLMM(object, null, ...)
\method{r.squaredGLMM}{merMod}(object, null, envir = parent.frame(), pj2014 = FALSE, ...)
}

\arguments{
  \item{object}{a fitted linear model object. }
  \item{null}{optionally, a null model, including only random effects. See \sQuote{Details}. }
  \item{envir}{optionally, the \code{environment} in which the null model is to be evaluated.
      Defaults to the current frame. See \code{\link{eval}}. }
  \item{pj2014}{logical, if \code{TRUE} and \code{object} is of \code{poisson} 
  family, the result will include \Rsqx{GLMM} using original formulation of 
  Johnson (2014). This requires fitting \code{object} with an observation-level
  random effect term added. }
  \item{\dots}{additional arguments, ignored.}
}

\value{
	\code{r.squaredGLMM} returns a two-column numeric \code{matrix}, each (possibly
	named) row holding values for marginal and conditional \Rsqx{GLMM}
	calculated with different methods, such as \dQuote{delta},
	\dQuote{log-normal}, \dQuote{trigamma}, or \dQuote{theoretical} for models
	of \code{binomial} family. 	
}

\details{

There are two types of \Rsqx{GLMM}: marginal and conditional.

\emph{Marginal \Rsqx{GLMM}} represents the variance explained by the fixed 
effects, and is defined as:

\mydequation{R_{GLMM(m)}^{2}= \frac{\sigma_f^2}{\sigma_f^2 + \sigma_{\alpha}^2 + \sigma_{\varepsilon }^2}
}{R_GLMM(m)² = (\sigma_f²) / (\sigma_f² + \sigma_\alpha² + \sigma_\epsilon²)
}{R_GLMM(m)^2 = (sigma_f^2) / (sigma_f^2 + sigma_alpha^2 + sigma_epsilon^2)
}

\emph{Conditional \Rsqx{GLMM}} represents the variance explained by the
entire model, including both fixed and random effects. It is calculated
by the equation:

\mydequation{R_{GLMM(c)}^{2}= \frac{\sigma_f^2 + \sigma_{\alpha}^2}{\sigma_f^2 + \sigma_{\alpha}^2 + \sigma_{\varepsilon }^2}
}{R_GLMM(c)² = (\sigma_f² + \sigma_\alpha²) / (\sigma_f² + \sigma_\alpha² + \sigma_\epsilon²)
}{R_GLMM(c)^2 = (sigma_f^2 + sigma_alpha^2) / (sigma_f^2 + sigma_alpha^2 + sigma_epsilon^2)
}

where \myequation{\sigma_f^2}{\sigma_f²}{sigma_f^2}
is the variance of the fixed effect components,
\myequation{\sigma_{\alpha}}{\sigma_\alpha²}{sigma_alpha^2} 
is the variance of the random effects, and  
\myequation{\sigma_\epsilon^2}{\sigma_\epsilon²}{sigma_epsilon^2}
is the \dQuote{observation-level} variance.

Three methods are available for deriving the observation-level variance
\eqn{\sigma_\varepsilon}: the delta method, lognormal approximation and using the
trigamma function.

The delta method can be used with for all distributions and link functions,
while lognormal approximation and trigamma function are limited to distributions
with logarithmic link. Trigamma-estimate is recommended whenever available.
Additionally, for binomial distributions, theoretical variances exist
specific for each link function distribution.

\emph{Null model}. Calculation of the observation-level variance involves in 
some cases fitting a \emph{null} model containing no fixed effects other than 
intercept, otherwise identical to the original model (including all the random 
effects). When using \code{r.squaredGLMM} for several models differing only in 
their fixed effects, in order to avoid redundant calculations, the null model 
object can be passed as the argument \code{null}.
Otherwise, a null model will be fitted \emph{via} updating the original model.
This assumes that all the variables used in the original model call have the 
same values as when the model was fitted. The function warns about this when 
fitting the null model is required. This warnings can be disabled by setting  
 \code{options(MuMIn.noUpdateWarning = TRUE)}. 

}

\note{

\strong{Important}: as of \pkg{MuMIn} version 1.41.0, 
\code{r.squaredGLMM} returns a revised statistics based on Nakagawa et 
al. (2017) paper. The returned value's format also has changed (it is a 
\code{matrix} rather than a numeric vector as before). Pre-1.41.0 version of the
function calculated the \dQuote{theoretical} \Rsqx{GLMM} for \code{binomial}
models.

\Rsqx{GLMM} can be calculated also for fixed-effect models. In
the simpliest case of \acronym{OLS} it reduces to \code{var(fitted) /
(var(fitted) + deviance / 2)}. Unlike likelihood-ratio based \Rsq for 
\acronym{OLS}, value of this statistic differs from that of
the classical \Rsq.

Currently methods exist for classes: \code{merMod}, \code{lme},
\code{glmmTMB}, \code{glmmADMB}, \code{glmmPQL}, \code{cpglm}(\code{m}) and
(\code{g})\code{lm}.

For families other than gaussian, Gamma, poisson, binomial and negative binomial,
the residual variance is obtained using \code{\link[insight]{get_variance}} 
from package \CRANpkg{insight}.

See note in \code{\link{r.squaredLR}} help page for comment on using \Rsq in 
model selection.

}

\author{
Kamil Barto\enc{ń}{n}. This implementation is based on \R code from 
\sQuote{Supporting Information} for Nakagawa et al. (2014), 
(the extension for random-slopes) Johnson (2014), and includes developments from
Nakagawa et al. (2017).
}

\references{
Nakagawa, S., Schielzeth, H. 2013 A general and simple method for obtaining
\Rsq from Generalized Linear Mixed-effects Models. \emph{Methods in
Ecology and Evolution} \bold{4}, 133--142.

Johnson, P. C. D. 2014 Extension of Nakagawa & Schielzeth’s \Rsqx{GLMM} to random 
slopes models. \emph{Methods in Ecology and Evolution} \bold{5}, 44--946.

Nakagawa, S., Johnson, P. C. D., Schielzeth, H. 2017 The coefficient of
determination \Rsq and intra-class correlation coefficient from generalized
linear mixed-effects models revisited and expanded.
\emph{J. R. Soc. Interface} \bold{14}, 20170213.

}

\seealso{
\code{\link{summary.lm}}, \code{\link{r.squaredLR}}

\code{\link[performance]{r2}} from package \CRANpkg{performance} calculates
\Rsqx{GLMM} also for variance at different levels, with optional confidence 
intervals. \CRANpkg{r2glmm} has functions for \Rsq and partial \Rsq.

}

\examples{

\dontshow{ if(require(nlme)) \{   }
data(Orthodont, package = "nlme")

fm1 <- lme(distance ~ Sex * age, ~ 1 | Subject, data = Orthodont)

fmnull <- lme(distance ~ 1, ~ 1 | Subject, data = Orthodont)

r.squaredGLMM(fm1)
r.squaredGLMM(fm1, fmnull)
r.squaredGLMM(update(fm1, . ~ Sex), fmnull)

r.squaredLR(fm1)
r.squaredLR(fm1, null.RE = TRUE)
r.squaredLR(fm1, fmnull) # same result

\dontrun{
if(require(MASS)) {
    fm <- glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID, 
        family = binomial, data = bacteria, verbose = FALSE)
    fmnull <- update(fm, . ~ 1)
    r.squaredGLMM(fm)

    # Include R2GLMM (delta method estimates) in a model selection table:
    # Note the use of a common null model
    dredge(fm, extra = list(R2 = function(x) r.squaredGLMM(x, fmnull)["delta", ]))
    
}
}
\dontshow{  \}  }
}

\keyword{models}
