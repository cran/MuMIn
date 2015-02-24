library(MuMIn)
library(mgcv)
options(na.action = "na.fail")
set.seed(2) ## simulate some data...

dat <- gamSim(1,n = 400, dist = "normal", scale = 2)

dat$fac <- gl(2, 200)
fm1 <- gam(y ~ s(x1, by = fac) + s(x2, by = fac), data = dat)

fm1$formula <- update.formula(fm1, . ~ . + s(x1) + s(x2))
getAllTerms(fm1)

AICc(fm1)

dredge(fm1, subset = !(`s(x1)` & `s(x1, by = fac)`) & !(`s(x2)` & `s(x2, by = fac)`) )



fm2 <- gam(y ~ x0, data = dat)
fm3 <- gam(y ~ 1, data = dat)



model.sel(fm1, fm2, fm3)

zd <- dredge(fm1)
(model.sel(zd))

