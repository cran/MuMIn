\name{MuMIn-package}
\alias{MuMIn-package}
\alias{MuMIn}
\docType{package}
\encoding{utf-8}
\title{Multi-model inference}
\description{
The package \pkg{MuMIn} contains functions to streamline information-theoretic
model selection and carry out model averaging based on information criteria.
}

\details{

The suite of functions includes:
\describe{
	\item{\lcode{dredge}}{performs automated model selection by generating
		subsets of the supplied \sQuote{global} model and optional
		choices of other model properties (such as different link functions).
		The set of models can be generated with \sQuote{all possible}
		combinations or tailored according to specified conditions.
		}
	\item{\lcode{model.sel}}{creates a model selection table from
		selected models. }
	\item{\lcode{model.avg}}{calculates model-averaged parameters,
		along with standard errors and confidence intervals.
        The \lxcode{predict}{=predict.averaging} method 
        produces model-averaged predictions. }
	\item{\lcode{AICc}}{calculates the second-order Akaike information
		criterion. Some other criteria are provided, see below. }
	\item{\lcode{stdize}, \lcode{stdizeFit}, \lcode{std.coef},
		\lcode{partial.sd}}{can be used to standardise data and
		model coefficients by standard deviation or partial standard deviation.
		}
}

For a complete list of functions, use \code{library(help = "MuMIn")}.

By default, \AICc is used to rank models and obtain model
weights, although any information criterion can be used. At least the 
following are currently implemented in \R:
\lcode{AIC} and \lcode{BIC} in package \pkg{stats}, and
\lcode{QAIC}, \lcode{QAICc}, \lcode{ICOMP},
\lcode{CAICF}, and \link{Mallows' Cp} in \pkg{MuMIn}. There is also a
\lcode{DIC} extractor for \acronym{MCMC} models and a \lcode{QIC} for
\acronym{GEE}.

Many model fitting functions in \R are supported. For a complete list,
see \link[=MuMIn-models]{the list of supported models}.

In addition to \dQuote{regular} information criteria, model averaging can be 
performed using various types of model weighting algorithms:
  \link[=BGWeights]{Bates-Granger},
  \link[=BGWeights]{bootstrapped},
  \link[=cos2Weights]{cos-squared},
  \link[=jackknifeWeights]{jackknife},
  \link[=stackingWeights]{stacking}, or
  \link[=armWeights]{ARM}. 
  These weighting functions are mainly applicable to \code{glm}s.

}

\author{ Kamil Barto\enc{ń}{n} }

\references{
Burnham, K. P. and Anderson, D. R. 2002 \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed. New York,
Springer-Verlag.
}


\seealso{
\lcode{AIC}, \lcode{step} or \lxcode{stepAIC}{MASS} for stepwise
model selection by \AIC.
}

\examples{
\dontshow{oop <- }
options(na.action = "na.fail") #  change the default "na.omit" to prevent models 
                               #  from being fitted to different datasets in 
                               #  case of missing values.

fm1 <- lm(y ~ ., data = Cement)
ms1 <- dredge(fm1)

# Visualize the model selection table:
\dontshow{ if(require(graphics)) \{ }
par(mar = c(3,5,6,4))
plot(ms1, labAsExpr = TRUE)
\dontshow{ \} }
model.avg(ms1, subset = delta < 4)

confset.95p <- get.models(ms1, cumsum(weight) <= .95)
avgmod.95p <- model.avg(confset.95p)
summary(avgmod.95p)
confint(avgmod.95p)
\dontshow{options(oop)}
}

\keyword{package}
\keyword{models}
