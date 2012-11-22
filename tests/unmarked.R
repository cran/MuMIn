if(length(.find.package("unmarked", quiet = TRUE)) != 0L) {

library(MuMIn)
library(stats4)
library(unmarked) # cheating R check
#do.call(library, alist(unmarked)) # cheating R check

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

MuMIn:::fixCoefNames(names(coef(fm2oc)))

#str(fm4oc)
#showMethods("backTransform")
#getMethod("backTransform", "unmarkedLinComb")
#backTransform(fm4oc)

# dredge(fm2oc, eval=F, fixed=~psi(habitat))

#(dd <- dredge(fm2oc, fixed = ~psi(habitat), trace = T))
(dd <- dredge(fm2oc, fixed = ~psi(habitat)))


model.sel(dd, rank = "AIC")
models <- get.models(dd[1:3])

summary(ma1 <- model.avg(models))
summary(ma2 <- model.avg(dd[1:3]))
summary(ma3 <- model.avg(model.sel(model.sel(dd, rank = "AIC"), rank = "AICc")[1:3]))

predict(ma1, type = "det")

stopifnot(!any(is.na(coefTable(ma1))))
stopifnot(!any(is.na(coefTable(ma2))))
stopifnot(!any(is.na(coefTable(ma3))))
stopifnot(isTRUE(all.equal(coefTable(ma1), coefTable(ma2))))
stopifnot(isTRUE(all.equal(coefTable(ma2), coefTable(ma3))))
stopifnot(!any(is.na(dd[, "psi(habitat)"])))

summary(model.avg(dd, delta <= 4))

#family.unmarkedFit <- function (object, ...) NA

# Model selection
print(mod.sel(fm1oc, fm2oc, fm3oc))

#models <- list(fm1oc, fm2oc, fm3oc)
#traceback()
#model.avg(fm1oc, fm2oc, fm3oc)
#model.avg(fm1oc, fm2oc)
#getAllTerms(fm1oc)

print(summary(model.avg(fm1oc, fm2oc, fm3oc)))

}
