library(MuMIn)
library(mgcv)
set.seed(2) ## simulate some data...
dat <- gamSim(1,n = 400, dist = "normal", scale = 2)
fm1 <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
fm2 <- gam(y ~ x0, data = dat)
fm3 <- gam(y ~ 1, data = dat)

model.sel(fm1, fm2, fm3)

zd <- dredge(fm1)
(model.sel(zd))

