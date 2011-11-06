# Test varia

require(MuMIn)

#print(packageDescription("MuMIn", fields = "Version"))


# TEST binary response ---------------------------------------------------------
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive=20-numdead)
budworm.lg <- glm(SF ~ sex*ldose, family=binomial)
dd <- dredge(budworm.lg, trace=FALSE)
gm <- get.models(dd, 1:4)
model.avg(gm)

# The same, but use cbind directly in the formula
budworm.lg <- glm(cbind(numdead, numalive=20-numdead) ~ sex*ldose, family=binomial)
dd <- dredge(budworm.lg, trace=TRUE)
gm <- get.models(dd, 1:4)
avgmod <- model.avg(gm)


# TEST for consistency of vcov and se calculation ------------------------------
if(!isTRUE(all.equal(avgmod$avg.model[,2],  sqrt(diag(vcov(avgmod))),  tolerance = .001)))
stop("'vcov' has a problem")



# TEST evaluation from within function -----------------------------------------

budworm <- data.frame(ldose = rep(0:5, 2), numdead = c(1, 4, 9, 12, 18, 20, 0,
	2, 6, 10, 12, 16), sex = factor(rep(c("M", "F"), c(6, 6))))
budworm$SF <- cbind(numdead = budworm$numdead, numalive = 20 - budworm$numdead)

# evaluate within an exotic environment
(function(dat) (function(dat2) {
	#mod <- glm(SF ~ sex*ldose, data = dat2, family = "quasibinomial", trace=T)
	mod <- glm(SF ~ sex*ldose, data = dat2, family = "binomial", trace=F)
	#mod <- glm(SF ~ sex*ldose, data = budworm, family = "binomial", trace=F)
	print(dd <- dredge(mod, rank = "QAIC", chat = summary(budworm.lg)$dispersion))
	gm <- get.models(dd, family="binomial")
	#print(sys.frames())
	summary(model.avg(gm))
})(dat))(budworm)


rm(list=ls())

### Bootstrap
##data(Cement)
##bCement <- Cement
##nobs <- nrow(Cement)
##gm1 <- lm(y~X1+X2+X3+X4, bCement)
##ms <- dredge(gm1, m.min=1)
##nrow(ms)
##cl <- attr(ms, "calls")
##bootCement2 <- t(sapply(1:5000, function(i) {
##	bCement <- Cement[sample(nobs, replace=TRUE), ]
##	models <- lapply(cl, eval, sys.frame(sys.nframe()))
##	m <- which.min(sapply(models, AICc))
##	matchCoef(models[[m]], gm1)
##}))
##bootCement <- rbind(bootCement2, bootCement)
##
##cbind(
##global = coef(gm1),
##mavg = coef(model.avg(get.models(ms, 1:5))),
##average = apply(bootCement, 2, mean, na.rm=T),
##se = apply(bootCement, 2, sd, na.rm=T) /
##	sqrt(colSums(!apply(bootCement, 2, is.na))),
##sd = apply(bootCement, 2, sd, na.rm=T)
##)
##model.avg(get.models(ms, 1:5))$avg.model


# END TESTS
