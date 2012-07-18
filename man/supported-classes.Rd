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
\item \code{gam}, \code{gamm} (\pkg{mgcv});
\item \code{gamm4} (\pkg{gamm4});
\item \code{glmmML} (\pkg{glmmML});
\item \code{glmmadmb} (\pkg{glmmADMB} from R-Forge);
\item \code{hurdle}, \code{zeroinfl} (\pkg{pscl});
\item \code{negbin}, \code{betabin} (class \code{glimML}, package \pkg{aod});
\item \code{sarlm}, \code{spautolm} (\pkg{spdep});
\item \code{coxph}, \code{survreg} (\pkg{survival});
\item \code{coxme}, \code{lmekin} (\pkg{coxme});
\item \code{rq} (\pkg{quantreg});

\item functions from package \pkg{unmarked}, (within the class \code{unmarkedFit})
}

\code{model.sel} and \code{dredge} accept also \code{MCMCglmm} models
	(package \pkg{MCMCglmm}) with \code{\link{DIC}} as the \code{rank} function.

Other classes are also likely to be supported, in particular if they inherit
from one of the above classes. In general, the models averaged with
\code{model.avg} may belong to different types, provided they use the same data
and response, and if it is valid to do so. This applies also to constructing
model selection tables with \code{model.sel}.

}

\seealso{
\code{\link{model.avg}}, \code{\link{model.sel}} and \code{\link{dredge}}.
}

\keyword{package}
