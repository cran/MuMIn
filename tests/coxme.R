library(coxme)
library(MuMIn)


temp <- with(lung, scale(cbind(age, wt.loss, meal.cal)))

rfit0 <- coxme(Surv(time, status) ~ ph.ecog * ph.karno + (age | 1) + (wt.loss | 1),
	data = lung)

stopifnot(is.numeric(coeffs(rfit0)))

dd <- dredge(rfit0, eval = T, trace = T)

model.sel(dd, rank = AIC)
summary(ma <- model.avg(dd))

#fit1 <- lme(effort ~ Type, data=ergoStool, random= ~1|Subject/ran1, method="ML")

if(exists("lmekin", mode = "function", envir = asNamespace("coxme"))) {
	fit2 <- lmekin(effort ~ Type + (1|Subject), data=ergoStool)
	dd <- dredge(fit2, trace = T)
	summary(ma <- model.avg(dd))
}
