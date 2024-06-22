\name{get.models}
\alias{get.models}
\alias{pget.models}
\encoding{utf-8}

\title{Retrieve models from selection table}

\description{
    Generate or extract a list of fitted model objects from a 
	\code{"model.selection"} table or component models from the averaged model
    (\code{"averaging"} object), optionally using parallel computation in a 
    cluster.
}


\usage{
get.models(object, subset, cluster = NA, ...)
}

\arguments{
    \item{object}{ object returned by \code{\link{dredge}}, 
        \code{\link{model.sel}} or \code{\link{model.avg}}.
    }
    \item{subset}{subset of models, an expression evaluated within the model
        selection table (see \sQuote{Details}).}
    \item{cluster}{ optionally, a \code{"cluster"} object. If it is a valid 
        cluster, models are evaluated using parallel computation. }
    \item{\dots}{additional arguments to update the models. For example, one 
        may want to fit models with \acronym{REML} (e.g. argument 
        \code{REML = TRUE} in some modelling functions) while using 
        \acronym{ML} for model selection.}
}

\value{
    \code{\link{list}} of fitted model objects.
}

\details{
The argument \code{subset} must be explicitely provided. This is to assure that 
a potentially long list of models is not fitted unintentionally. To evaluate all
 models, set \code{subset} to \code{NA} or \code{TRUE}. 

If \code{subset} is a character vector, it is interpreted as names of rows to be
selected.

}

\note{

\code{"model.selection"} tables created by \code{model.sel} or averaged models
 created by \code{model.avg} from a list of model objects (as opposed to those
 created with model selection tables) store the component models as part of the
 object - in these cases \code{get.models}  simply extracts the items from 
 these lists. Otherwise the models have to be fitted. Therefore, using 
 \code{get.models} following \code{dredge} is not efficient as the 
 requested models are fitted twice. If the number of generated models is 
 reasonable, consider using \code{lapply(dredge(..., evaluate = FALSE), eval)}, 
 which generates a list of all model calls and evaluates them into a list of 
 model objects. 

Alternatively, \code{getCall} and \code{eval} can be used to compute a model out of the 
\code{"model.selection"} table (e.g. \code{eval(getCall(<model.selection>, i))}, where 
\code{i} is the model index or name). 

\code{pget.models} is still available, but is deprecated.
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
\code{\link{dredge}} and \code{\link{pdredge}}, \code{\link{model.avg}}

\code{\link{makeCluster}} in packages \pkg{parallel} and \pkg{snow}

}

\examples{
# Mixed models:

\dontshow{ if(require(nlme)) \{ }
fm2 <- lme(distance ~ age + Sex, data = Orthodont,
    random = ~ 1 | Subject, method = "ML")
ms2 <- dredge(fm2)

# Get top-most models, but fitted by REML:
(confset.d4 <- get.models(ms2, subset = delta < 4, method = "REML"))

\dontrun{
# Get the top model:
get.models(ms2, subset = 1)[[1]]
}

\dontshow{  \} }
}

\keyword{models}
