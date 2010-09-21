# Test support for different classes of models

require(MuMIn)

# TEST gls --------------------------------------------------------------------------------
library(nlme)
data(Ovary, package = "nlme")
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
		   correlation = corAR1(form = ~ 1 | Mare))
dd <- dredge(fm1, trace=T)
gm <- get.models(dd, 1:4)
ma <- model.avg(gm)
predict(ma)
predict(ma, data.frame(Mare=1, Time=range(Ovary$Time)))

detach(package:nlme)
rm(list=ls())

# TEST nlme -------------------------------------------------------------------------------
library(nlme)
data(Orthodont, package = "nlme")

fm2 <- lme(distance ~ Sex*age + age*Sex, data = Orthodont,
		   random = ~ 1|Subject / Sex,
		   method = "ML")
dd <- dredge(fm2, trace=T)
gm <- get.models(dd, 1:4)
ma <- model.avg(gm)
predict(ma)
predict(ma, data.frame(Sex="Male", Subject="M01", age=8:12))

detach(package:nlme)
rm(list=ls())

# TEST lmer -------------------------------------------------------------------------------
library(lme4)
data(Orthodont, package = "nlme")

fm2 <-lmer(distance ~ Sex*age + (1|Subject) + (1|Sex), data = Orthodont)
dd <- dredge(fm2, trace=T)
gm <- get.models(dd, 1:4)
ma <- model.avg(gm)
#predict(ma)
#predict(ma, data.frame(Sex="Male", Subject="M01", age=8:12))

detach(package:lme4)
rm(list=ls())

# TEST lm ---------------------------------------------------------------------------------
data(Orthodont, package = "nlme")

fm1 <- lm(distance ~ Sex*age + age*Sex, data = Orthodont)
dd <- dredge(fm1, trace=T)
gm <- get.models(dd, 1:4)
ma <- model.avg(gm)
predict(ma)
predict(ma, data.frame(Sex="Male", age=8:12))

rm(list=ls())

# TEST glm --------------------------------------------------------------------------------
data(Cement, package = "MuMIn")

nseq <- function(x, len=length(x)) seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE),
	length=len)

fm1 <- glm(y ~ (X+X1+X2+X3)^3, data = Cement)
dd <- dredge(fm1, trace=T)
gm <- get.models(dd, 1:10)
ma <- model.avg(gm)
predict(ma) == predict(ma, Cement)
predict(ma, lapply(Cement, nseq))

rm(list=ls())

# TEST rlm --------------------------------------------------------------------------------
library(MASS)
data(Cement, package = "MuMIn")

nseq <- function(x, len=length(x)) seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE),
	length=len)

fm1 <- rlm(y ~ X+X1+X2*X3, data = Cement)
dd <- dredge(fm1, trace=T)
gm <- get.models(dd, 1:10)
ma <- model.avg(gm)
predict(ma) == predict(ma, Cement)
predict(ma, lapply(Cement, nseq))

rm(list=ls())
detach(package:MASS)

# TEST multinom --------------------------------------------------------------------------------
library(nnet)
library(MASS)

example(birthwt)
bwt.mu <- multinom(low ~ ., bwt)

dd <- dredge(bwt.mu, trace=T)
gm <- get.models(dd, 1:10)
ma <- model.avg(gm)
#predict(bwt.mu)
#predict(ma) // Cannot average factors!

rm(list=ls())
detach(package:nnet)

# TEST gam --------------------------------------------------------------------------------
require(mgcv)

dat <- gamSim(1, n = 500, dist="poisson", scale=0.1)

gam1 <- gam(y ~ s(x0) + s(x1) + s(x2) +  s(x3) + (x1+x2+x3)^2,
	family = poisson, data = dat, method = "REML")


dd <- dredge(gam1, subset=(!`s(x1)` | !x1) & (!`s(x2)` | !x2) & (!`s(x3)` | !x3),
			 fixed="x1", trace=T)
gm <- get.models(dd, cumsum(weight) <= .95)
ma <- model.avg(gm)
predict(ma)

rm(list=ls())
detach(package:mgcv)

# TEST spautolm ---------------------------------------------------------------------------
library(spdep)
example(NY_data)

esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="full", verbose=FALSE)

dd <- dredge(esar1f)
gm <- get.models(dd, cumsum(weight) <= .98)
ma <- model.avg(gm)
resid(ma)


rm(list=ls())
detach(package:spdep)

# TEST spautolm ---------------------------------------------------------------------------
library(spdep)
data(oldcol)

COL.errW.eig <- errorsarlm(CRIME ~ INC * HOVAL * OPEN, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen", quiet=TRUE)

dd <- dredge(COL.errW.eig)
gm <- get.models(dd, cumsum(weight) <= .98)
ma <- model.avg(gm)

predict(ma)
formula(ma)
resid(ma)

rm(list=ls())
detach(package:spdep)

# END TESTS