###
# Example of model selection with models from 'unmarked' package
# with parallel execution
###

require(MuMIn)
library(unmarked)
require(stats4)
require(parallel) || require(snow)

# Set up the cluster
ncores <- if(exists("detectCores", mode = "function")) 
	detectCores() else getOption("cl.cores", 2)
clust <- makeCluster(ncores)

data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site, 
    obsCovs = mallard.obs)
	
# Fit the global model
(ufm.mallard <- pcount(~ ivel + date + I(date^2) ~ length + elev + forest, 
    mallardUMF, K = 30))

invisible(clusterEvalQ(clust, library(unmarked, logical = TRUE)))
# 'stats4' is needed for AIC to work with unmarkedFit objects but is not loaded
# automatically with 'unmarked'.
invisible(clusterEvalQ(clust, library(stats4, logical = TRUE)))

clusterExport(clust, "mallardUMF")
invisible(clusterCall(clust, "library", "stats4", character.only = TRUE))

# For comparison, single-threaded run:
#system.time(print(pdd1 <- pdredge(ufm.mallard,
#   subset = `p(date)` | !`p(I(date^2))`, rank = AIC)))

system.time(pdd2 <- pdredge(ufm.mallard, clust,
    subset = `p(date)` | !`p(I(date^2))`, rank = AIC, extra = "adjR^2"))


# select the top models and null model
subset(pdd2, delta < 2 | df == min(df))

# Remove the warnings permanently
attr(pdd2, "warnings") <- NULL

# Compare with the model selection table from unmarked
# the statistics should be identical:
models <- pget.models(pdd2, clust, delta < 2 | df == min(df))

modSel(fitList(fits = structure(models, names = model.names(models,
    labels = getAllTerms(ufm.mallard)))), nullmod = "(Null)")
	
stopCluster(clust)

########################