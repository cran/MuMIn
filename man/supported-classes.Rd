\name{MuMIn-models}
\alias{MuMIn-models}
\encoding{utf-8}
\title{List of supported models}
\description{
List of model classes accepted by \code{model.avg}, \code{model.sel},
	and \code{dredge}.
}

\details{
Fitted model objects that can be used with model selection and model averaging
functions include those returned by:

\itemize{
\item \code{lm}, \code{glm} (package \pkg{stats});
\item \code{rlm}, \code{glm.nb} and \code{polr} (\pkg{MASS});
\item \code{multinom} (\pkg{nnet});
\item \code{lme}, \code{gls} (\pkg{nlme});
\item \code{lmer}, \code{glmer} (\pkg{lme4});
\item \code{cpglm}, \code{cpglmm} (\pkg{cplm});
\item \code{gam}, \code{gamm}* (\pkg{mgcv});
\item \code{gamm4}* (\pkg{gamm4});
\item \code{glmmML} (\pkg{glmmML});
\item \code{glmmadmb} (\pkg{glmmADMB} from R-Forge);
\item \code{asreml} (\pkg{asreml}, non-free commercial package; REML comparisons only);
\item \code{hurdle}, \code{zeroinfl} (\pkg{pscl});
\item \code{negbin}, \code{betabin} (class \code{glimML}, package \pkg{aod});
\item \code{aodml}, \code{aodql} (\pkg{aods3});
\item \code{betareg} (\pkg{betareg});
\item \code{sarlm}, \code{spautolm} (\pkg{spdep});
\item \code{spml}* (if fitted by \acronym{ML}, \pkg{splm});
\item \code{coxph}, \code{survreg} (\pkg{survival});
\item \code{coxme}, \code{lmekin} (\pkg{coxme});
\item \code{rq} (\pkg{quantreg});
\item \code{clm} and \code{clmm} (\pkg{ordinal});
\item \code{logistf} (\pkg{logistf}); %% Bias-reduced logistic regression
\item \code{crunch}*, \code{pgls} (\pkg{caper});
\item functions from package \pkg{unmarked} (within the class \code{unmarkedFit});
\item \code{mark} and related functions (class \code{mark} from package
	\pkg{RMark}). Note currently \code{dredge} can only manipulate \code{formula}
	element of the argument \code{model.parameters}, keeping its other elements
	intact.
}


Generalized Estimation Equation model implementations: \code{\link[geepack]{geeglm}} 
from package \pkg{geepack}, \code{\link[gee]{gee}} from \pkg{gee}, and 
\code{yags} from \pkg{yags} (from R-Forge) can be used with \code{\link{QIC}} 
as the selection criterion.

\code{MCMCglmm}* models (package \pkg{MCMCglmm}) with e.g. \code{\link{DIC}} as the 
\code{rank} function are accepted by \code{model.sel} and \code{dredge}.

Other classes are also likely to be supported, in particular if they inherit
from one of the above classes. In general, the models averaged with
\code{model.avg} may belong to different types (e.g. \code{glm} and \code{gam}), 
provided they use the same data and response, and if it is valid to do so. 
This applies also to constructing model selection tables with \code{model.sel}.
}

\note{
* In order to use \code{gamm}, \code{gamm4}, \code{spml (> 1.0.0)},
  \code{crunch} or \code{MCMCglmm} with \code{dredge}, an
  \code{\link{updateable}} wrapper for these functions should be created.
}

\seealso{
\code{\link{model.avg}}, \code{\link{model.sel}} and \code{\link{dredge}}.
}

\keyword{package}
