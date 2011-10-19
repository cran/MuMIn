library(MuMIn)

# rm(list=ls(all=TRUE))

.checkPkg <- function(package) length(.find.package(package, quiet=TRUE)) != 0
if(.checkPkg("unmarked")) {
library(unmarked)

# Simulate occupancy data
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
for(i in 1:nSites) {
    y[i,] <- rbinom(nReps, 1, Z[i]*p[i])
    }

umf <- unmarkedFrameOccu(y = y, siteCovs = covariates)
#head(umf)
#summary(umf)


# Fit some models
fm1oc <- occu(~1 ~1, umf)
fm2oc <- occu(~veght+habitat ~veght+habitat, umf)
fm3oc <- occu(~veght ~veght+habitat, umf)


# Model selection
print(mod.sel(fm1oc, fm2oc, fm3oc))
print(summary(model.avg(fm1oc, fm2oc, fm3oc)))

}
