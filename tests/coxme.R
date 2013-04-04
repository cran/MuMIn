if(length(find.package("coxme", quiet = TRUE)) == 2) {

library(coxme)
library(MuMIn)


lung$temp <- with(lung, scale(cbind(age, wt.loss, meal.cal)))
lung1 <- na.omit(lung)

rfit0 <- coxme(Surv(time, status) ~ ph.ecog * ph.karno + (age | 1) + (wt.loss | 1),
	data = lung1)

stopifnot(is.numeric(coeffs(rfit0)))

dd <- dredge(rfit0, eval = T, trace = T)
coeffs(dd)
model.sel(dd, rank = AIC)
summary(ma <- model.avg(dd))

library(survival)

fm0 <- coxph(formula = Surv(time, status) ~ 1, data = lung1)
fm1 <- coxph(formula = Surv(time, status) ~ temp, data = lung1)
fme <- coxme(formula = Surv(time, status) ~ ph.ecog + (temp | 1), data = lung1)
fme1 <- coxme(formula = Surv(time, status) ~ ph.ecog + (age | 1) + (wt.loss | 1), data = lung1)
fme2 <- coxme(formula = Surv(time, status) ~ ph.ecog + (wt.loss | 1), data = lung1)

model.sel(fm0, fm1, fme, fme1, fme2)


#fit1 <- lme(effort ~ Type, data=ergoStool, random= ~1|Subject/ran1, method="ML")

if(exists("lmekin", mode = "function", envir = asNamespace("coxme"))) {
	fit2 <- lmekin(effort ~ Type + (1|Subject), data=ergoStool)
	dd <- dredge(fit2, trace = T)
	summary(ma <- model.avg(dd))
}

}