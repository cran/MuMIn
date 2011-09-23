# Test support for different classes of models

require(MuMIn)
.checkPkg <- function(package) length(.find.package(package, quiet=TRUE)) != 0

# TEST gls --------------------------------------------------------------------------------
if(.checkPkg("nlme")) {
library(nlme)
data(Ovary, package = "nlme")
fm1 <- gls(follicles ~ sin(2*pi*Time) + cos(2*pi*Time), Ovary,
		   correlation = corAR1(form = ~ 1 | Mare))
dd <- dredge(fm1, trace=T)
gm <- get.models(dd, 1:4)
ma <- model.avg(gm, revised=T)

dd
mod.sel(gm)

summary(ma)
confint(ma)

#ma <- model.avg(gm, revised=F)

predict(ma)
predict(ma, data.frame(Mare=1, Time=range(Ovary$Time)))

detach(package:nlme); rm(list=ls())
}

#dredge(fm1, rank=BIC)
#dredge(fm1, rank=AIC)

# TEST (quasi)poisson --------------------------------------------------------------------

d.AD <- data.frame(counts =c(18,17,15,20,10,20,25,13,12), outcome = gl(3,1,9),
treatment = gl(3,3))
glm.qD93 <- glm(counts ~ outcome + treatment, family=quasipoisson(), data=d.AD)
glm.D93 <- glm(counts ~ outcome + treatment, family=poisson(), data=d.AD)

dd <- dredge(glm.qD93)
summary(model.avg(dd, subset= delta <= 10))
dd <- dredge(glm.D93)
summary(model.avg(dd, subset= delta <= 10))

subset(dd, delta <= 10)
mod.sel(get.models(dd, subset=delta <= 10))


rm(list=ls())

# TEST glmmML --------------------------------------------------------------------

if(.checkPkg("glmmML")) {

library(glmmML)
dat <- data.frame(y = rbinom(100, prob = rep(runif(20), rep(5, 20)), size = 1),
	x = rnorm(100), x2 = rnorm(100), id = factor(rep(1:20, rep(5, 20))))

fm1 <- glmmML(y ~ x*x2, data = dat, cluster = id)
summary(model.avg(dredge(fm1), subset = delta <= 4))

detach(package:glmmML); rm(list=ls())

}

# TEST nlme --------------------------------------------------------------------
if(.checkPkg("nlme")) {
library(nlme)
data(Orthodont, package = "nlme")

#:: Model-averaging mixed models :::::::::::::::::::::::::::::::::::::::::::::::
# Fitting by REML
fm2 <- lme(distance ~ Sex*age + age*Sex, data = Orthodont,
		   random = ~ 1|Subject / Sex, method = "REML")

# Model selection: ranking by AICc which uses ML
dd <- dredge(fm2, trace=T, rank="AICc", REML=FALSE)

#attr(dd, "rank.call")

# Get models (which are fitted by REML, like the global model)
gm <- get.models(dd, 1:4)

#mod.sel(gm)

##::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#maML <- model.avg(gm, method = "NA", rank="AICc", rank.args = list(REML=FALSE))
#maREML <- model.avg(gm, method = "NA", rank="AICc", rank.args = list(REML=TRUE))

#summary(maML)
#summary(maREML)

#ma <- model.avg(gm, revised = F)
summary(maNA <- model.avg(gm, revised = T, method = "NA"))
summary(ma0 <- model.avg(gm, revised = T, method = "0"))
confint(maNA)
confint(ma0)

#dredge(fm2, rank=BIC)
predict(ma0, data.frame(Sex="Male", Subject="M01", age=8:12))

detach(package:nlme); rm(list=ls())
}

# TEST lmer -------------------------------------------------------------------------------
if(.checkPkg("lme4")) {

library(lme4)
data(Orthodont, package = "nlme")

fm2 <- lmer(distance ~ Sex*age + (1|Subject) + (1|Sex), data = Orthodont)

dd <- dredge(fm2, trace=T)
gm <- get.models(dd, 1:4)
(ma <- model.avg(gm))


#predict(ma)
#predict(ma, data.frame(Sex="Male", Subject="M01", age=8:12))

# TEST dredge with update'd model
dd <- dredge(update(fm2, REML=FALSE), trace=T)

# update.mer does not expand dots, so here we have a call:
# lmer(formula = distance ~ Sex + (1 | Subject), data = Orthodont,
#    REML = ..2, model = ..3)
dd <- dredge(update(fm2, REML=F, model=F), trace=T)

detach(package:lme4); rm(list=ls())
}

# TEST lm ---------------------------------------------------------------------------------
if(.checkPkg("nlme")) {
data(Orthodont, package = "nlme")

fm1 <- lm(distance ~ Sex*age + age*Sex, data = Orthodont)
dd <- dredge(fm1, trace=T)
gm <- get.models(dd, 1:4)
ma <- model.avg(gm, revised=F)

vcov(ma)
summary(ma)
confint(ma)

predict(ma)
predict(ma, se.fit=TRUE)
predict(ma, data.frame(Sex="Male", age=8:12))

rm(list=ls())
}

# TEST glm --------------------------------------------------------------------------------
data(Cement, package = "MuMIn")

#invisible(lapply(dir("e:/Dokumenty/R-Forge/mumin/pkg/R", pattern="\\.R$", full.names = T), source, keep.source =F, echo=F))

nseq <- function(x, len=length(x)) seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE),
	length=len)

fm1 <- glm(y ~ (X+X1+X2+X3)^2, data = Cement)
dd <- dredge(fm1, trace=T)
#dd <- dredge(fm1)

gm <- get.models(dd, 1:10)
#ma <- model.avg(gm)
ma <- model.avg(gm, method="NA")
ma <- model.avg(gm, method="0")
vcov(ma)

summary(ma)
predict(ma) == predict(ma, Cement)
predict(ma, se.fit=T)
predict(ma, lapply(Cement, nseq))


rm(list=ls())
# TEST rlm --------------------------------------------------------------------------------

if (.checkPkg("MASS")) {
library(MASS)
data(Cement, package = "MuMIn")

nseq <- function(x, len=length(x)) seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE),
	length=len)

fm1 <- rlm(y ~ X+X1+X2*X3, data = Cement)
dd <- dredge(fm1, trace=T)
gm <- get.models(dd, 1:10)
ma <- model.avg(gm)
stopifnot(all(predict(ma) == predict(ma, Cement)))
predict(ma, lapply(Cement, nseq, len=30), se.fit=TRUE)

rm(list=ls()); detach(package:MASS)

# TEST multinom --------------------------------------------------------------------------------
if (.checkPkg("nnet")) {
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
summary(ma)
#predict(bwt.mu)
# predict(ma) // Cannot average factors!

rm(list=ls()); detach(package:nnet)
}}


# TEST gam --------------------------------------------------------------------------------
if (.checkPkg("mgcv")) {

suppressPackageStartupMessages(library(mgcv))
RNGkind("Mersenne")
set.seed(0) ## simulate some data...
dat <- gamSim(1,n=400,dist="normal",scale=2)
#gam1 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat)

ops <- options(warn = -1)

gam1 <- gam(y ~ s(x0) + s(x1) + s(x2) +  s(x3) + (x1+x2+x3)^2,
	data = dat, method = "ML")

dd <- dredge(gam1, subset=!`s(x0)` & (!`s(x1)` | !x1) & (!`s(x2)` | !x2) & (!`s(x3)` | !x3), fixed="x1", trace=T)

gm <- get.models(dd, cumsum(weight) <= .95)
ma <- model.avg(gm)

predict(ma, dat[1:10, ], se.fit=T)
options(ops)

rm(list=ls()); detach(package:mgcv)
}

# TEST spautolm ---------------------------------------------------------------------------

if (.checkPkg("spdep"))
if(!is.null(tryCatch(suppressPackageStartupMessages(library(spdep)), error=function(e) NULL))) {

suppressMessages(example(NY_data, echo = FALSE))

esar1f <- spautolm(Z ~ PEXPOSURE * PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="full", verbose=FALSE)

dd <- dredge(esar1f, m.max=3, fixed=~PEXPOSURE)
gm <- get.models(dd, cumsum(weight) <= .99)
ma <- model.avg(gm)
summary(ma)
signif(resid(ma), 5)[1:10]

rm(list=ls())

# TEST spautolm ---------------------------------------------------------------------------
#suppressPackageStartupMessages(library(spdep))
data(oldcol)

COL.errW.eig <- errorsarlm(CRIME ~ INC * HOVAL * OPEN, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen", quiet=TRUE)

dd <- dredge(COL.errW.eig)
gm <- get.models(dd, cumsum(weight) <= .98)
ma <- model.avg(gm)

# XXX: Error in t(vapply(models, function(m) {: error in evaluating the argument 'x' in selecting a method for function 't': Error in m.tTable[, 2L] : incorrect number of dimensions
summary(ma)
predict(ma)[1:10]

rm(list=ls()); detach(package:spdep)

} # library(spdep)



# TEST glm.nb ---------------------------------------------------------------------------
if (.checkPkg("MASS")) {
require(MASS)

quine.nb1 <- glm.nb(Days ~ 0+Sex/(Age + Eth*Lrn), data = quine)
#quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)

dredge(quine.nb1) # Wrong
dredge(quine.nb1, marg.ex="Sex") # Right
ma <- model.avg(dredge(quine.nb1, marg.ex="Sex"), subset=cumsum(weight)<=.9999)

pred <- predict(ma, se=TRUE)
#pred <- cbind(pred$fit, pred$fit - (2 * pred$se.fit), pred$fit + (2 * pred$se.fit))
#matplot(pred, type="l")
#matplot(family(quine.nb1)$linkinv(pred), type="l")

rm(list=ls()); detach(package:MASS)
}

# TEST quasibinomial ---------------------------------------------------------------------------

budworm <- data.frame(ldose = rep(0:5, 2), numdead = c(1, 4, 9, 13, 18, 20, 0,
	2, 6, 10, 12, 16), sex = factor(rep(c("M", "F"), c(6, 6))))
budworm$SF = cbind(numdead = budworm$numdead, numalive = 20 - budworm$numdead)

budworm.lg <- glm(SF ~ sex*ldose + sex*I(ldose^2), data = budworm, family = quasibinomial)
budworm.lg1 <- glm(SF ~ sex*ldose + sex*I(ldose^2), data = budworm, family = binomial)

dd <- dredge(budworm.lg, rank = "QAIC", chat = summary(budworm.lg)$dispersion)
dd <- dredge(budworm.lg) # should be the same
mod <- get.models(dd, seq(nrow(dd)))

# Note: this works:
# ma <- model.avg(mod)
# but this will not: ('rank' attribute passed from 'dredge' is lost)
# ma <- model.avg(mod)
# so, need to supply them
ma <- model.avg(mod[1:5], rank="QAICc", rank.args = list(chat = 0.403111))

rm(list=ls());

# TEST polr {MASS} ---------------------------------------------------------------------------
#if (.checkPkg("MASS")) {
#library(MASS)
#library(MuMIn)
#options(contrasts = c("contr.treatment", "contr.poly"))
#house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)

#dd <- dredge(house.plr)
#mod <- get.models(dd, 1:6)
#model.avg(mod)
#}

# END TESTS
