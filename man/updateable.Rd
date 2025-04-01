\name{updateable}
\alias{updateable}
\alias{updateable2}
\alias{get_call}
\alias{uGamm}
\alias{fitdistr2}
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

## updateable wrapper for MASS::fitdistr
fitdistr2(x, densfun, start, ...)

}

\arguments{
    \item{FUN}{function to be modified, found \emph{via} \lcode{match.fun}. }
	\item{eval.args}{optionally a character vector of function arguments' names
	      to be evaluated in the stored call. See \sQuote{Details}. }
	
	\item{Class}{optional character vector naming class(es) to be set onto the
		result of \code{FUN} (not possible if the result is an S4 object). }
	
	\item{x}{for \code{get_call}, an object from which the call should be 
        extracted. For \code{fitdistr2}, a numeric vector passed to 
        \lxcode{fitdistr}{MASS}. }
	
	\item{formula, random}{arguments to be passed to \code{gamm} or \code{gamm4}}
	\item{lme4}{if \code{TRUE}, \code{gamm4} is called, \code{gamm} otherwise. }
    
    \item{densfun, start}{Arguments passed to \lxcode{fitdistr}{MASS}. }

    \item{\dots}{Arguments passed to respective wrapped functions. }

}

\details{

	Most model fitting functions in \R return an object that can be updated or
	re-fitted \emph{via} \lcode{update}. This is possible thanks to the 
    function call stored in the object, which can be used (possibly modified)
    later on. It is also used by \code{dredge} to generate submodels. Some
    functions (such as \code{mgcv::gamm} or \code{MCMCglmm::MCMCglmm}) do not 
    provide their result with the \code{call} element. To work around this, 
    \code{updateable} can be used on such a function to store the call. The 
    resulting \dQuote{wrapper} should be used in exactly the same way as the 
    original function.
    
    \code{updateable} can also be used to repair an existing \code{call} element,
	e.g. if it contains \link[=dots]{dotted names} that prevent re-evaluation
	of a call.
    
    The \code{eval.args} argument specifies the names of the function arguments 
    to be evaluated in the stored call. This is useful if, for example, the 
    model object does not have a \code{formula} element or does not 
    store the formula in any other way, and the modelling function has been 
    called with the formula specified as the variable name. In this case, the 
    default \code{formula} method will try to retrieve the formula from the 
    stored \code{call}, which does not guarantee that the variable will be 
    available at the time of retrieval, or that the value of that variable will
    be the same as that used to fit the model (this is demonstrated in the 
    last \sQuote{example}).

}

\value{
   \code{updateable} returns a function with the same arguments as \code{FUN}, 
   wrapping a call to  \code{FUN} and adding an element named \code{call} to its
   result if possible, otherwise an attribute \code{"call"} (if the returned 
   value is atomic or an S4 object).
}

\note{
  \code{get_call} is similar to \lcode{getCall} (defined in package 
  \pkg{stats}), but it can also extract the \code{call} when it is an 
  \lxcode{attribute}{=attr} (and not an element of the object). Because the
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
\lcode{update}, \lcode{getCall}, \lcode{getElement},
\lcode{attributes}

\lxcode{gamm}{mgcv}, \lxcode{gamm4}{gamm4}

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
