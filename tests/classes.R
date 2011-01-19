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

detach(package:nlme); rm(list=ls())

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

detach(package:nlme); rm(list=ls())

# TEST lmer -------------------------------------------------------------------------------

library(lme4)
data(Orthodont, package = "nlme")

fm2 <- lmer(distance ~ Sex*age + (1|Subject) + (1|Sex), data = Orthodont)

dd <- dredge(fm2, trace=T)
gm <- get.models(dd, 1:4)
ma <- model.avg(gm)
#predict(ma)
#predict(ma, data.frame(Sex="Male", Subject="M01", age=8:12))

# TEST dredge with update'd model
dd <- dredge(update(fm2, REML=FALSE), trace=T)

# update.mer does not expand dots, so here we have a call:
# lmer(formula = distance ~ Sex + (1 | Subject), data = Orthodont,
#    REML = ..2, model = ..3)
dd <- dredge(update(fm2, REML=F, model=F), trace=T)

detach(package:lme4); rm(list=ls())

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

fm1 <- glm(y ~ (X+X1+X2+X3)^2, data = Cement)
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

rm(list=ls()); detach(package:MASS)

# TEST multinom --------------------------------------------------------------------------------
library(nnet); library(MASS)

# Trimmed-down model from example(birthwt)
data(birthwt)

bwt <- with(birthwt, data.frame(
		low = low,
		race = factor(race, labels = c("white", "black", "other")),
		ptd = factor(ptl > 0),
		smoke = (smoke > 0)
		))

options(contrasts = c("contr.treatment", "contr.poly"))
bwt.mu <- multinom(low ~ ., data = bwt)

dd <- dredge(bwt.mu, trace=T)
gm <- get.models(dd, seq(nrow(dd)))
ma <- model.avg(gm)
#predict(bwt.mu)
# predict(ma) // Cannot average factors!

rm(list=ls()); detach(package:nnet)

# TEST gam --------------------------------------------------------------------------------
suppressPackageStartupMessages(library(mgcv))

RNGkind("Mersenne")
set.seed(8)
dat <- gamSim(1, n = 50, dist="poisson", scale=0.1)
gam1 <- gam(y ~ s(x0) + s(x1) + s(x2) +  s(x3) + (x1+x2+x3)^2,
	family = poisson, data = dat, method = "REML")
dd <- dredge(gam1, subset=!`s(x0)` & (!`s(x1)` | !x1) & (!`s(x2)` | !x2) & (!`s(x3)` | !x3),
			 fixed="x1", trace=T)
gm <- get.models(dd, cumsum(weight) <= .95)
ma <- model.avg(gm)
predict(ma)

rm(list=ls()); detach(package:mgcv)

# TEST spautolm ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(spdep))
suppressMessages(example(NY_data, echo = FALSE))

esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="full", verbose=FALSE)

dd <- dredge(esar1f)
gm <- get.models(dd, cumsum(weight) <= .98)
ma <- model.avg(gm)
signif(resid(ma), 5)

rm(list=ls()); detach(package:spdep)

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
signif(resid(ma), 5)

rm(list=ls()); detach(package:spdep)

# TEST glm.nb ---------------------------------------------------------------------------
require(MASS)

quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)

dredge(quine.nb1) # Wrong
dd <- dredge(quine.nb1, marg.ex="Sex") # Right
model.avg(get.models(dd, 1:5))

rm(list=ls()); detach(package:MASS)

# TEST quasibinomial ---------------------------------------------------------------------------

budworm <- data.frame(ldose = rep(0:5, 2), numdead = c(1, 4, 9, 13, 18, 20, 0,
	2, 6, 10, 12, 16), sex = factor(rep(c("M", "F"), c(6, 6))))
budworm$SF = cbind(numdead = budworm$numdead, numalive = 20 - budworm$numdead)

budworm.lg <- glm(SF ~ sex*ldose + sex*I(ldose^2), data = budworm, family = quasibinomial)

dd <- dredge(budworm.lg, rank = "QAIC", chat = summary(budworm.lg)$dispersion)
dd <- dredge(budworm.lg) # should be the same
mod <- get.models(dd, seq(nrow(dd)))

# Note: this works:
# ma <- model.avg(mod)
# but this will not: ('rank' attribute passed from 'dredge' is lost)
# ma <- model.avg(mod)
# so, need to supply them
ma <- model.avg(mod[1:5], rank="QAICc", rank.args = list(chat = 0.403111))

rm(list=ls())

# TEST polr {MASS} ---------------------------------------------------------------------------

#library(MASS)
#library(MuMIn)
#options(contrasts = c("contr.treatment", "contr.poly"))
#house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)

#dd <- dredge(house.plr)
#mod <- get.models(dd, 1:6)
#model.avg(mod)


# END TESTS
