%% \newcommand{\lcode}{\code{\link{#1}}}


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
\item \lxcode{rlm}{MASS}, \lxcode{glm.nb}{MASS} and \lxcode{polr}{MASS} (\CRANpkg{MASS});
\item \lxcode{multinom}{nnet} (\CRANpkg{nnet});
\item \lxcode{lme}{nlme}, \lxcode{gls}{nlme} (\CRANpkg{nlme});
\item \lxcode{lmer}{lme4}, \lxcode{glmer}{lme4} (\CRANpkg{lme4});
\item \lxcode{cpglm}{cplm}, \lxcode{cpglmm}{cplm} (\CRANpkg{cplm});
\item \lxcode{gam}{mgcv}, \lxcode{gamm}{mgcv}* (\CRANpkg{mgcv});
\item \lxcode{gamm4}{gamm4}* (\CRANpkg{gamm4});
\item \lxcode{gamlss}{gamlss} (\CRANpkg{gamlss});
\item \lxcode{glmmML}{glmmML} (\CRANpkg{glmmML});
\item \code{glmmadmb} (\href{http://glmmadmb.r-forge.r-project.org/}{\pkg{glmmADMB}}
	from R-Forge);
\item \lxcode{glmmTMB}{glmmTMB} (\CRANpkg{glmmTMB});
\item \lxcode{MCMCglmm}{MCMCglmm}* (\CRANpkg{MCMCglmm});
\item \code{asreml} (non-free commercial package \pkg{asreml}; allows only for 
	\acronym{REML} comparisons);
\item \lxcode{hurdle}{pscl}, \lxcode{zeroinfl}{pscl} (\CRANpkg{pscl});
\item \lxcode{negbin}{aod}, \lxcode{betabin}{aod} (class \code{"glimML"}), 
	package \CRANpkg{aod});
\item \unlxcode{aodml}{aods3}, \unlxcode{aodql}{aods3} (\unCRANpkg{aods3});
\item \lxcode{betareg}{betareg} (\CRANpkg{betareg});
\item \lxcode{brglm}{brglm} (\CRANpkg{brglm});
\item \code{*sarlm} models, \lxcode{spautolm}{spatialreg} (\CRANpkg{spatialreg});
\item \lxcode{spml}{splm}* (if fitted by \acronym{ML}, \CRANpkg{splm});
\item \lxcode{coxph}{survival}, \lxcode{survreg}{survival} (\CRANpkg{survival});
\item \lxcode{coxme}{coxme}, \lxcode{lmekin}{coxme} (\CRANpkg{coxme});
\item \lxcode{rq}{quantreg} (\CRANpkg{quantreg});
\item \lxcode{clm}{ordinal} and \lxcode{clmm}{ordinal} (\CRANpkg{ordinal});
\item \lxcode{logistf}{logistf} (\CRANpkg{logistf}); %% Bias-reduced logistic regression
\item \lxcode{crunch}{caper}*, \lxcode{pgls}{caper} (\CRANpkg{caper});
\item \lxcode{maxlike}{maxlike} (\CRANpkg{maxlike});
\item most \code{"unmarkedFit"} objects from package \CRANpkg{unmarked});
\item \lxcode{mark}{RMark} and related functions (class \code{mark} from package
	\CRANpkg{RMark}). Currently \code{dredge} can only manipulate \code{formula}
	element of the argument \code{model.parameters}, keeping its other elements
	intact;
\item \lxcode{fitdistr}{MASS} mostly useful for model selection with 
    \code{model.sel}. Use of \lxcode{fitdistr2}{=updateable} wrapper 
    function is recommended.
}

Generalized Estimation Equation model implementations: 
\lxcode{geeglm}{geepack} from package \CRANpkg{geepack}, 
\lxcode{gee}{gee} from \CRANpkg{gee}, 
\lxcode{geem}{geeM} from \CRANpkg{geeM}, 
\unlxcode{wgee}{wgeesel} from \unCRANpkg{wgeesel}, 
and \code{yags} from \href{http://yags.r-forge.r-project.org/}{\pkg{yags}} (on 
R-Forge) can be used with \lcode{QIC} as the selection criterion.


Further classes may also be supported, in particular if they inherit from one of
the classes listed above. In general, models averaged using \code{model.avg} can
belong to different types (e.g. \code{glm} and \code{gam}), provided they use 
the same data and response, and, obviously, if it is valid to do so. This also 
applies to the construction of model selection tables using \code{model.sel}.

}

\note{
* In order to use \code{gamm}, \code{gamm4}, \code{spml (> 1.0.0)},
  \code{crunch} or \code{MCMCglmm} with \code{dredge}, an
  \lcode{updateable} wrapper for these functions should be created.
}

\seealso{
\lcode{model.avg}, \lcode{model.sel} and \lcode{dredge}.
}

\keyword{package}
