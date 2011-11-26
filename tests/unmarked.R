library(MuMIn)
library(unmarked)

# Simulate occupancy data
set.seed(1)
umfOccu <- local({
	nSites <- 100
	nReps <- 5
	covariates <- data.frame(veght=rnorm(nSites),
		habitat=factor(c(rep('A', 50), rep('B', 50))))
	psipars <- c(-1, 1, -1)
	ppars <- c(1, -1, 0)
	X <- model.matrix(~veght+habitat, covariates) # design matrix
	psi <- plogis(X %*% psipars)
	p <- plogis(X %*% ppars)
	y <- matrix(NA, nSites, nReps)
	Z <- rbinom(nSites, 1, psi)       # true occupancy state
	for(i in 1:nSites) y[i,] <- rbinom(nReps, 1, Z[i]*p[i])
	unmarkedFrameOccu(y = y, siteCovs = covariates)
})

# Fit some models
fm1oc <- occu(~1 ~1, umfOccu)
fm2oc <- occu(~veght+habitat ~veght*habitat, umfOccu)
fm3oc <- occu(~habitat ~veght+habitat, umfOccu)
fm4oc <- occu(~veght ~veght+habitat, umfOccu)

# dredge(fm2oc, eval=F, fixed=~psi(habitat))

(dd <- dredge(fm2oc, fixed=~psi(habitat)))

#attr(dd, "terms")
#coef(fm4oc)

stopifnot(!any(is.na(dd[, "psi(habitat)"])))

summary(model.avg(dd, delta <= 4))

# Model selection
print(mod.sel(fm1oc, fm2oc, fm3oc))
print(summary(model.avg(fm1oc, fm2oc, fm3oc)))
