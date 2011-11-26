# Test support for different classes of models

require(MuMIn)
.checkPkg <- function(package) length(.find.package(package, quiet=TRUE)) != 0

# TEST gls --------------------------------------------------------------------------------
#if(.checkPkg("nlme")) {
library(nlme)

fm1Dial.gls <- gls(rate ~(pressure + I(pressure^2) + I(pressure^3) + I(pressure^4))*QB,
      Dialyzer, method="ML")

varying <- list(
	correlation = alist(
		AR1_0.771=corAR1(0.771, form = ~ 1 | Subject), 	AR1=corAR1(), NULL),
	weights = alist(vp.press=varPower(form = ~ pressure), NULL)
	)

dd <- dredge(fm1Dial.gls, m.max = 1, m.min = 1, fixed=~pressure, varying =
	varying, extra = "R^2")
models <- get.models(dd, 1:4)
ma <- model.avg(models, revised=T)

summary(ma)
predict(ma)
predict(ma, data.frame(QB = 300, pressure = seq(0.24, 3, length = 10)))

# example(corGaus)
fm1BW.lme <- lme(weight ~ Time * Diet, data = BodyWeight, random = ~ Time, method="ML")

varying <- list(
	correlation = alist(corExp(form = ~ Time), corGaus(form = ~ Time)),
	weights = alist(NULL, varPower())
)

dd <- dredge(fm1BW.lme, m.max=1, fixed=~Time, varying=varying)
#dd <- dredge(fm1, trace=T)
models <- get.models(dd, 1:4)
ma <- model.avg(models, revised=T)
mod.sel(models)
summary(ma)
confint(ma)
#ma <- model.avg(gm, revised=F)

detach(package:nlme); rm(list=ls())
#}

#dredge(fm1, rank=BIC)
#dredge(fm1, rank=AIC)

# TEST (quasi)poisson --------------------------------------------------------------------

d.AD <- data.frame(counts =c(18,17,15,20,10,20,25,13,12), outcome = gl(3,1,9),
treatment = gl(3,3))
#fm2 <- glm(counts ~ outcome + treatment, family=quasipoisson(), data=d.AD)
fm2 <- glm(counts ~ outcome + treatment, family=poisson(), data=d.AD)

#dd <- dredge(fm2)
#summary(model.avg(dd, subset= delta <= 10))
dd <- dredge(fm2, rank=QAIC, chat=deviance(fm2) / df.residual(fm2))
dd <- dredge(fm2)
summary(model.avg(dd, subset= delta <= 10))

dredge(fm2, rank=QAIC, chat=deviance(fm2) / df.residual(fm2))

subset(dd, delta <= 10)
mod.sel(get.models(dd, subset = delta <= 10))


rm(list=ls())

# TEST glmmML --------------------------------------------------------------------

if(.checkPkg("glmmML")) {
#require(MuMIn)
library(glmmML)
set.seed(100)
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
#maML <- model.avg(gm, rank="AICc", rank.args = list(REML=FALSE))
#maREML <- model.avg(gm, rank="AICc", rank.args = list(REML=TRUE))
#summary(maML)
#summary(maREML)

#ma <- model.avg(gm, revised = F)
summary(ma <- model.avg(gm, revised = T))
confint(ma)

#dredge(fm2, rank=BIC)
predict(ma, data.frame(Sex="Male", Subject="M01", age=8:12))

detach(package:nlme); rm(list=ls())
}

# TEST lmer -------------------------------------------------------------------------------
if(.checkPkg("lme4")) {

set.seed(1)
library(lme4)
data(Orthodont, package = "nlme")

Orthodont$rand <- runif(nrow(Orthodont))
fm2 <- lmer(log(distance) ~ rand*Sex*age + (1|Subject) + (1|Sex), data = Orthodont, REML=FALSE)

#fm0 <- lmer(distance ~ 1 + (1|Subject) + (1|Sex), data = Orthodont, REML=FALSE)
#fm00 <- lm(distance ~ 1 , data = Orthodont)

dd <- dredge(fm2, trace=T)
gm <- get.models(dd, 1:4)
(ma <- model.avg(gm))

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

dispersion <- function(object) {
	wts <- weights(object); if(is.null(wts)) wts <- 1
	sum((wts * resid(object, type="working")^2)[wts > 0]) / df.residual(object)
}


dd <- dredge(fm1, extra = alist(dispersion))
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

nseq <- function(x, len=length(x)) seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE),
	length=len)

fm1 <- glm(y ~ (X1+X2+X3+X4)^2, data = Cement)
dd <- dredge(fm1, trace=T)

gm <- get.models(dd, 1:10)

ma <- model.avg(gm)
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

fm1 <- rlm(y ~X1+X2*X3+X4, data = Cement)
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

options(warn=1)
dd <- dredge(esar1f, m.max=1,  fixed=~PEXPOSURE,
	varying = list(
		family = list("CAR", "SAR"),
		method=list("Matrix_J", "full")
	), trace=FALSE)
options(warn=0)


#dd <- dredge(esar1f, m.max=3, fixed=~PEXPOSURE)
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

# XXX: Error in t(vapply(models, function(m) {: error in evaluating the argument
# 'x' in selecting a method for function 't': Error in m.tTable[, 2L] :
# incorrect number of dimensions
summary(ma)
predict(ma)[1:10]

rm(list=ls()); detach(package:spdep)

} # library(spdep)

# TEST glm.nb ---------------------------------------------------------------------------
if (.checkPkg("MASS")) {
require(MASS)

quine.nb1 <- glm.nb(Days ~ 0 + Sex/(Age + Eth*Lrn), data = quine)
#quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)

ms <- dredge(quine.nb1, marg.ex = "Sex")
models <- get.models(ms, 1:5)
summary(model.avg(models))

dredge(quine.nb1) # Wrong
dredge(quine.nb1, marg.ex="Sex") # Right
ma <- model.avg(dredge(quine.nb1, marg.ex="Sex"), subset=cumsum(weight)<=.9999)

pred <- predict(ma, se=TRUE)
#pred <- cbind(pred$fit, pred$fit - (2 * pred$se.fit), pred$fit + (2 * pred$se.fit))
#matplot(pred, type="l")
#matplot(family(quine.nb1)$linkinv(pred), type="l")

rm(list=ls()); detach(package:MASS)
}

# TEST quasibinomial -----------------------------------------------------------

budworm <- data.frame(ldose = rep(0:5, 2), numdead = c(1, 4, 9, 13, 18, 20, 0,
	2, 6, 10, 12, 16), sex = factor(rep(c("M", "F"), c(6, 6))))
budworm$SF = cbind(numdead = budworm$numdead, numalive = 20 - budworm$numdead)

qbinomial <- function(...) {
	res <- quasibinomial(...)
	res$aic <- binomial(...)$aic
	res
}
budworm.qqlg <- glm(SF ~ sex*ldose + sex*I(ldose^2), data = budworm, family = qbinomial)
#budworm.qlg <- glm(SF ~ sex*ldose + sex*I(ldose^2), data = budworm, family = quasibinomial)
budworm.lg <- glm(SF ~ sex*ldose + sex*I(ldose^2), data = budworm, family = binomial)

dd <- dredge(budworm.lg, rank = "QAIC", chat = summary(budworm.lg)$dispersion)
#dd <- dredge(budworm.lg) # should be the same
mod <- get.models(dd, seq(nrow(dd)))

# Note: this works:
# ma <- model.avg(mod)
# but this will not: ('rank' attribute passed from 'dredge' is lost)
# ma <- model.avg(mod)
# so, need to supply them
ma <- model.avg(mod[1:5], rank="QAICc", rank.args = list(chat = 0.403111))

rm(list=ls())

# TEST polr {MASS} -------------------------------------------------------------
#if (.checkPkg("MASS")) {
#library(MASS)
#library(MuMIn)
#options(contrasts = c("contr.treatment", "contr.poly"))
#house.plr <- polr(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)

#dd <- dredge(house.plr)
#mod <- get.models(dd, 1:6)
#model.avg(mod)
#}


# TEST coxph -------------------------------------------------------------------

library(survival)

bladder1 <- bladder[bladder$enum < 5, ]

fmcph <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) +  cluster(id), bladder1)

r.squared.coxph <- function(object, ...) {
	logtest <- -2 * (object$loglik[1L] - object$loglik[2L])
	c(rsq = 1 - exp(-logtest/object$n), maxrsq = 1 - exp(2 * object$loglik[1L]/object$n))
}

ms <- dredge(fmcph, fixed=c("cluster(id)", "strata(enum)"), extra = "r.squared.coxph")
fits <- get.models(ms, delta < 5)
summary(model.avg(fits))

fmsrvrg <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='weibull',
    scale=1)
summary(model.avg(dredge(fmsrvrg), delta  < 4))
# nobs(fmsrvrg)

rm(list=ls())
detach(package:survival)

# END TESTS

###model1 <- glm(cbind(ncases, ncontrols) ~ agegp + tobgp * alcgp,
###              data = esoph, family = binomial())
###model1 <- glm(cbind(ncases, ncontrols) ~ agegp,
###              data = esoph, family = binomial())
###model2 <- update(model1, contrasts = list(agegp="contr.treatment"))
###model3 <- update(model1, .~. + tobgp, contrasts = list(agegp="contr.treatment"))
###model3$contrasts
###
###models <- list(model1, model2, model3)
###x <- rbindDataFrameList(lapply(sapply(models, getElement, "contrasts"), as.data.frame))
###
###as.list(x)
###
###any(apply(x, 2, function(x) length(na.omit(unique(x)))) > 1)
###
####system.time(for(i in 1:5000) lapply(lapply(lapply(x, unique), na.omit), length))
####system.time(for(i in 1:5000) lapply(x, function(x) length(na.omit(unique(x)))))
###
###summary(model.avg(model1, model2, model3))
###
###
###model1$contrasts
###
###model2$contrasts
###
###
###ms1 <- dredge(model1, rank=Cp)
###
####plot(dredge(model1, rank=ICOMP))
###dredge(model1, rank=BIC, extra="R^2")
###
###model1$contrasts
###summary(model.avg(ms1, delta <=10))
###
###summary(model2 <- glm(case ~ age+parity+education+spontaneous+induced,
###                data=infert, family = binomial()))
###
###dredge(model2, rank=BIC, extra="nobs")

# require(survival)
# model3 <- clogit(case~spontaneous+induced+strata(stratum),data=infert)
# model3 <- clogit(case~spontaneous+strata(stratum),data=infert)

#require(quantreg)
#data(stackloss)
#fm1 <- rq(stack.loss ~ Water.Temp + Acid.Conc. + Air.Flow, tau = .5, data = stackloss)
#ms <- dredge(fm1, varying = list(tau=list(.75, .5, .25)), extra="R^2")
#fm1 <- lm(stack.loss ~ Water.Temp + Acid.Conc. + Air.Flow,  data = stackloss)
#models <- get.models(ms, delta <= 4)
#model.avg(ms, delta <= 4)
#plot(ms)

#
#library(brglm)
#data(lizards)
#lizg <- glm(cbind(grahami, opalinus) ~ (height + diameter +
#    light + time)^2, family = binomial(logit), data=lizards)
#
#glm(cbind(grahami, opalinus) ~ (height + diameter +
#    light + time)^3, family = poisson("log"), data=lizards)
#dd <- dredge(lizg)
#plot(dd, col=6:7)

#brglm(cbind(grahami, opalinus) ~ height + diameter + light + time,
#	family = binomial(logit), data=lizards, method = "brglm.fit")
#
#lizards.brglm <- brglm(cbind(grahami, opalinus) ~ height + diameter +
#                  light + time, family = binomial(logit), data=lizards,
#                  method = "brglm.fit")
