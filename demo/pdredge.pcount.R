###
# Example of model selection with models from 'unmarked' package
# with parallel execution
###
require(parallel)
library(MuMIn)
library(unmarked)

# Set up the cluster
ncores <- if(exists("detectCores", mode = "function"))
	detectCores() else getOption("cl.cores", 2)
clust <- try(makeCluster(getOption("cl.cores", 2), type = "PSOCK"))
if(!inherits(clust, "cluster")) stop("Could not set up the cluster")

data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
    obsCovs = mallard.obs)

# Fit the global model
(ufm.mallard <- pcount(~ ivel + date + I(date^2) ~ length + elev + forest,
    mallardUMF, K = 30))

invisible(clusterEvalQ(clust, library(unmarked, logical = TRUE)))

clusterExport(clust, "mallardUMF")

# For comparison, single-threaded run:
#system.time(print(pdd1 <- pdredge(ufm.mallard,
#   subset = `p(date)` | !`p(I(date^2))`, rank = AIC)))

system.time(pdd2 <-
	dredge(ufm.mallard, cluster = clust,
    subset = (`p(date)` || !`p(I(date^2))`),
	rank = AIC, extra = "adjR^2", eval = TRUE))


# select the top models and null model
subset(pdd2, delta < 2 | df == min(df))

# Remove the warnings permanently
attr(pdd2, "warnings") <- NULL

# Compare with the model selection table from 'unmarked'.
# The statistics should be identical:
models <- get.models(pdd2, delta < 2 | df == min(df), cluster = clust)

modSel(fitList(fits = structure(models, names = model.names(models,
    labels = getAllTerms(ufm.mallard)))), nullmod = "(Null)")

stopCluster(clust)
########################
