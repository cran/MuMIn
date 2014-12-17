\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}} 
\newcommand{\lxcode}{\code{\link[#2]{#1}}}
\newcommand{\lcode}{\code{\link{#1}}}


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
functions include those produced by:

\itemize{
\item \lcode{lm}, \lcode{glm} (package \pkg{stats});
\item \lcode{rlm}, \lcode{glm.nb} and \lcode{polr} (\CRANpkg{MASS});
\item \lcode{multinom} (\CRANpkg{nnet});
\item \lxcode{lme}{nlme}, \lxcode{gls}{nlme} (\CRANpkg{nlme});
\item \lxcode{lmer}{lme4}, \lxcode{glmer}{lme4} (\CRANpkg{lme4});
\item \lxcode{cpglm}{cplm}, \lxcode{cpglmm}{cplm} (\CRANpkg{cplm});
\item \lxcode{gam}{mgcv}, \lxcode{gamm}{mgcv}* (\CRANpkg{mgcv});
\item \lxcode{gamm4}{gamm4}* (\CRANpkg{gamm4});
\item \lxcode{glmmML}{glmmML} (\CRANpkg{glmmML});
\item \code{glmmadmb} (\href{http://glmmadmb.r-forge.r-project.org/}{\pkg{glmmADMB}}
	from R-Forge);
\item \code{asreml} (non-free commercial package \pkg{asreml}; allows only for 
	\acronym{REML} comparisons);
\item \lxcode{hurdle}{pscl}, \lxcode{zeroinfl}{pscl} (\CRANpkg{pscl});
\item \lxcode{negbin}{aod}, \lxcode{betabin}{aod} (class \code{"glimML"}), 
	package \CRANpkg{aod});
\item \lxcode{aodml}{aods3}, \lxcode{aodql}{aods3} (\CRANpkg{aods3});
\item \lxcode{betareg}{betareg} (\CRANpkg{betareg});
\item \lxcode{brglm}{brglm} (\CRANpkg{brglm});
\item \code{*sarlm} models, \lxcode{spautolm}{spdep} (\CRANpkg{spdep});
\item \lxcode{spml}{splm}* (if fitted by \acronym{ML}, \CRANpkg{splm});
\item \lxcode{coxph}{survival}, \lxcode{survreg}{survival} (\CRANpkg{survival});
\item \lxcode{coxme}{coxme}, \lxcode{lmekin}{coxme} (\CRANpkg{coxme});
\item \lxcode{rq}{quantreg} (\CRANpkg{quantreg});
\item \lxcode{clm}{ordinal} and \lxcode{clmm}{ordinal} (\CRANpkg{ordinal});
\item \lxcode{logistf}{logistf} (\CRANpkg{logistf}); %% Bias-reduced logistic regression
\item \lxcode{crunch}{caper}*, \lxcode{pgls}{caper} (\CRANpkg{caper});
\item \lxcode{maxlike}{maxlike} (\CRANpkg{maxlike});
\item functions from package \CRANpkg{unmarked} (within the class 
	\code{"unmarkedFit"});
\item \lxcode{mark}{RMark} and related functions (class \code{mark} from package
	\CRANpkg{RMark}). Currently \code{dredge} can only manipulate \code{formula}
	element of the argument \code{model.parameters}, keeping its other elements
	intact.
}

Generalized Estimation Equation model implementations: \lxcode{geeglm}{geepack}
from package \CRANpkg{geepack}, \lxcode{gee}{gee} from \CRANpkg{gee}, and 
\code{yags} from \href{http://yags.r-forge.r-project.org/}{\pkg{yags}} (from 
R-Forge) can be used with \code{\link{QIC}} as the selection criterion.

\lxcode{MCMCglmm}{MCMCglmm}* models (package \CRANpkg{MCMCglmm}) with e.g. 
\code{\link{DIC}} as the \code{rank} function are accepted by \code{model.sel} 
and \code{dredge}.

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
