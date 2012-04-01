\name{pdredge}
\alias{pdredge}

\encoding{utf-8}
\title{Automated model selection using parallel computation}
\description{
Parallelized version of \code{dredge}.
}

\usage{
pdredge(global.model, cluster = NA, beta = FALSE, evaluate = TRUE,
    rank = "AICc", fixed = NULL, m.max = NA, m.min = 0, subset,
    marg.ex = NULL, trace = FALSE, varying, extra, check = FALSE,
    ...)
}

\arguments{
    \item{global.model, beta, evaluate, rank}{
        see \code{\link{dredge}}. }
    \item{fixed, m.max, m.min, subset, marg.ex, varying, extra, ...}{
        see \code{\link{dredge}}. }
    \item{trace}{ displays the generated calls, but may not work as expected
        since the models are evaluated in batches rather than one by one. }
    \item{cluster}{ either a valid \code{cluster} object, or \code{NA} for a
        single threaded execution. }
    \item{check}{ logical, whether to evaluate the \code{global.model} in the
        cluster and compare with the original one using \code{all.equal}. }

}

\details{
All the dependencies for fitting the \code{global.model}, including the data,
	and any objects the modelling function will use, must be exported
	(\emph{via} e.g. \code{clusterExport}) into the clurter worker nodes, as
	well as the required packages must be loaded thereinto (e.g. \emph{via}
	\code{clusterEvalQ(..., library(package))}, before the cluster is used by
	\code{pdredge}.

	\code{pdredge} tries to check whether all the variables and functions used
	in the call to \code{global.model} are present in the cluster nodes'
	\code{.GlobalEnv} before proceeding further, and also, if \code{check} is
	\code{TRUE}, it will compare the \code{global.model} updated at the cluster
	nodes with the one given as argument. This additional test will, however,
    slow down the execution.

Use of \code{pdredge} should be considered mainly with large datasets and
	complex models, for which \code{dredge} takes a long time to complete.
	Otherwise there may be no perceptible improvement or even the parallel
	version may perform worse than a single-threaded one.
}


\value{
 See \code{\link{dredge}}.
}

\author{Kamil Barto\enc{ń}{n}}

\seealso{
    \code{makeCluster} and other cluster related functions in packages
    \pkg{parallel} or \pkg{snow}.
}


\examples{

\dontshow{
# Normally this should be simply "require(parallel) || require(snow)",
# but here we resort to an (ugly) trick to avoid MuMIn's dependency on one of
# these packages and still pass R-check:
if(MuMIn:::.parallelPkgCheck(quiet = TRUE)) \{
}

# One of these packages is required:
\dontrun{require(parallel) || require(snow)}

# From example(Beetle)
data(Beetle)

Beetle100 <- Beetle[sample(nrow(Beetle), 100, replace = TRUE),]

fm1 <- glm(Prop ~ dose + I(dose^2) + log(dose) + I(log(dose)^2),
    data = Beetle100, family = binomial)

msubset <- expression(xor(dose, `log(dose)`) & (dose | !`I(dose^2)`)
    & (`log(dose)` | !`I(log(dose)^2)`))
varying.link <- list(family = alist(logit = binomial("logit"),
    probit = binomial("probit"), cloglog = binomial("cloglog") ))

# Set up the cluster
clusterType <- if(length(.find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 2), type = clusterType))
\dontshow{if(inherits(clust, "cluster")) \{ }
clusterExport(clust, "Beetle100")

# noticeable gain only when data has about 3000 rows (Windows 2-core machine)
print(system.time(dredge(fm1, subset = msubset, varying = varying.link)))
print(system.time(pdredge(fm1, cluster = FALSE, subset = msubset,
    varying = varying.link)))
print(system.time(pdd <- pdredge(fm1, cluster = clust, subset = msubset,
    varying = varying.link)))

print(pdd)

\dontrun{
# A time consuming example with 'unmarked' model, based on example(pcount).
# Having enough patience, you can run this with 'demo(pdredge.pcount)'.
library(unmarked)
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
    obsCovs = mallard.obs)
(ufm.mallard <- pcount(~ ivel + date + I(date^2) ~ length + elev + forest,
    mallardUMF, K = 30))
clusterEvalQ(clust, library(unmarked))
clusterExport(clust, "mallardUMF")

# 'stats4' is needed for AIC to work with unmarkedFit objects but is not
# loaded automatically with 'unmarked'.
require(stats4)
invisible(clusterCall(clust, "library", "stats4", character.only = TRUE))

#system.time(print(pdd1 <- pdredge(ufm.mallard,
#   subset = `p(date)` | !`p(I(date^2))`, rank = AIC)))

system.time(print(pdd2 <- pdredge(ufm.mallard, clust,
    subset = `p(date)` | !`p(I(date^2))`, rank = AIC, extra = "adjR^2")))


# best models and null model
subset(pdd2, delta < 2 | df == min(df))

# Compare with the model selection table from unmarked
# the statistics should be identical:
models <- pget.models(pdd2, clust, delta < 2 | df == min(df))

modSel(fitList(fits = structure(models, names = model.names(models,
    labels = getAllTerms(ufm.mallard)))), nullmod = "(Null)")
}

stopCluster(clust)
\dontshow{
\} else # if(! inherits(clust, "cluster"))
message("Could not set up the cluster")
\}
}

}

\keyword{models}
