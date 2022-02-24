\name{updateable}
\alias{updateable}
\alias{updateable2}
\alias{get_call}
\alias{uGamm}
\alias{MuMIn-gamm}
\alias{gamm-wrapper}
\encoding{utf-8}

\title{Make a function return updateable result}
\description{
Creates a function wrapper that stores a call in the object returned by its
argument \code{FUN}.
}

\usage{
updateable(FUN, eval.args = NULL, Class)

get_call(x)

## updateable wrapper for mgcv::gamm and gamm4::gamm4
uGamm(formula, random = NULL, ..., lme4 = inherits(random, "formula"))
}

\arguments{
    \item{FUN}{function to be modified, found \emph{via} \code{\link{match.fun}}. }
	\item{eval.args}{optionally a character vector of function arguments' names
	      to be evaluated in the stored call. See \sQuote{Details}. }
	
	\item{Class}{optional character vector naming class(es) to be set onto the
		result of \code{FUN} (not possible with formal S4 objects). }
	
	\item{x}{an object from which the call should be extracted. }
	
	\item{formula, random, \dots}{arguments to be passed to \code{gamm} or \code{gamm4}}
	\item{lme4}{if \code{TRUE}, \code{gamm4} is called, \code{gamm} otherwise. }
}

\details{
	Most model fitting functions in \R return an object that can be updated or
	re-fitted \emph{via} \code{\link{update}}. This is thanks to the \code{call}
	stored in the object, which can be used (possibly modified) later on. It is
	also utilised by \code{dredge} to generate sub-models. Some functions (such
	as \code{gamm} or \code{MCMCglmm}) do not provide their result with the
	\code{call} element. To work that around, \code{updateable} can be used on
	that function to store the call. The resulting \dQuote{wrapper} should be
	used in exactly the same way as the original function.
	
	\code{updateable} can also be used to repair an existing \code{call} element,
	e.g. if it contains \link[=dots]{dotted names} that prevent re-evaluation
	of such a call.
	
	Argument \code{eval.args} specifies names of function arguments that should
	be evaluated in the stored call. This is useful when, for example, the model
	object does not have \code{formula} element. The default \code{formula} 
	method tries to retrieve formula from the stored \code{call},
	which works unless the formula has been given as a variable and value of
	that variable changed since the model was fitted (the last \sQuote{example}
	demonstrates this).
	
	
}

\value{
   \code{updateable} returns a function with the same arguments as \code{FUN}, 
   wrapping a call to  \code{FUN} and adding an element named \code{call} to its
   result if possible, otherwise an attribute \code{"call"} (if the returned 
   value is atomic or a formal S4 object).
}

\note{
  \code{get_call} is similar to \code{\link{getCall}} (defined in package 
  \pkg{stats}), but it can also extract the \code{call} when it is an 
  \code{\link[=attr]{attribute}} (and not an element of the object). Because the
  default \code{getCall} method cannot do that, the default \code{update} method
  will not work with atomic or S4 objects resulting from \code{updateable} 
  wrappers.

 \code{uGamm} sets also an appropriate class onto the result (\code{"gamm4"}
 and/or \code{"gamm"}), which is needed for some generics defined in \pkg{MuMIn}
 to work (note that unlike the functions created by \code{updateable} it has no
 formal arguments of the original function). As of version 1.9.2,
 \code{MuMIn::gamm} is no longer available.

}


\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{update}}, \code{\link{getCall}}, \code{\link{getElement}},
\code{\link{attributes}}

\code{\link[mgcv]{gamm}}, \code{\link[gamm4]{gamm4}}

}

\examples{

# Simple example with cor.test:

# From example(cor.test)
x <- c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
y <- c( 2.6,  3.1,  2.5,  5.0,  3.6,  4.0,  5.2,  2.8,  3.8)

ct1 <- cor.test(x, y, method = "kendall", alternative = "greater")

uCor.test <- updateable(cor.test)

ct2 <- uCor.test(x, y, method = "kendall", alternative = "greater")

getCall(ct1) # --> NULL
getCall(ct2)

#update(ct1, method = "pearson") --> Error
update(ct2, method = "pearson")
update(ct2, alternative = "two.sided")


## predefined wrapper for 'gamm':
\dontshow{ if(require(mgcv)) \{ }
set.seed(0)
dat <- gamSim(6, n = 100, scale = 5, dist = "normal")

fmm1 <- uGamm(y ~s(x0)+ s(x3) + s(x2), family = gaussian, data = dat, 
    random = list(fac = ~1))

getCall(fmm1)
class(fmm1)
\dontshow{ \} }
###

\dontrun{
library(caper)
data(shorebird)
shorebird <- comparative.data(shorebird.tree, shorebird.data, Species)

fm1 <- crunch(Egg.Mass ~ F.Mass * M.Mass, data = shorebird)

uCrunch <- updateable(crunch)

fm2 <- uCrunch(Egg.Mass ~ F.Mass * M.Mass, data = shorebird)

getCall(fm1)
getCall(fm2)
update(fm2) # Error with 'fm1'
dredge(fm2)
}

###
\dontrun{
# "lmekin" does not store "formula" element 
library(coxme)
uLmekin <- updateable(lmekin, eval.args = "formula")

f <- effort ~ Type + (1|Subject)
fm1 <- lmekin(f, data = ergoStool)
fm2 <- uLmekin(f, data = ergoStool)

f <- wrong ~ formula # reassigning "f"

getCall(fm1) # formula is "f"
getCall(fm2) 

formula(fm1) # returns the current value of "f" 
formula(fm2)
}
}

\keyword{utils}
