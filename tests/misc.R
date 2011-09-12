# Test varia

require(MuMIn)

print(packageDescription("MuMIn", fields = "Version"))


# TEST binary response ---------------------------------------------------------
ldose <- rep(0:5, 2)
numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
sex <- factor(rep(c("M", "F"), c(6, 6)))
SF <- cbind(numdead, numalive=20-numdead)
budworm.lg <- glm(SF ~ sex*ldose, family=binomial)
dd <- dredge(budworm.lg, trace=TRUE)
gm <- get.models(dd, 1:4)
model.avg(gm)

# The same, but use cbind directly in the formula
budworm.lg <- glm(cbind(numdead, numalive=20-numdead) ~ sex*ldose, family=binomial)
dd <- dredge(budworm.lg, trace=TRUE)
gm <- get.models(dd, 1:4)
model.avg(gm)
# TEST evaluation from within function -----------------------------------------

budworm <- data.frame(ldose = rep(0:5, 2), numdead = c(1, 4, 9, 13, 18, 20, 0,
	2, 6, 10, 12, 16), sex = factor(rep(c("M", "F"), c(6, 6))))
budworm$SF <- cbind(numdead = budworm$numdead, numalive = 20 - budworm$numdead)


(function(dat) (function(dat2) {
	mod <- glm(SF ~ sex*ldose, data = budworm, family = "quasibinomial", trace=T)
	mod1 <- glm(SF ~ sex*ldose, data = budworm, family = "binomial", trace=T)
	dd <- dredge(mod, rank = "QAIC", chat = summary(budworm.lg)$dispersion)
	gm <- get.models(dd, family="binomial")
	#print(sys.frames())
	model.avg(gm)
})(dat))(budworm)

# dredge(mod1, rank="QAIC", chat=summary(mod1)$dispersion)
# dredge(mod1, rank="qAIC", dispersion=summary(mod1)$dispersion)

# 1:100

# Sys.sleep(1)
# ll <- logLik(mod)

# ll2 <- sum(mod$y * log(mod$fitted.values) - mod$fitted.values - lgamma(mod$y+1))
# ll <- do.call("structure", c(ll2, attributes(ll)))


# bbmle::qAIC(mod, dispersion=summary(mod)$dispersion)


# summary(mod1)$dispersion
# fm1

# for(i in 1:10) log(-1)
# logLik(fm1)

# ?deviance
# (fm1)

# traceback()

# logLik(mod1)

# deviance(mod)
# deviance(mod1)

# glm.fit



# QAIC(mod, chat=)
# AIC(ll)

# AIC(mod1)

# logLik(mod1)

# # predictors being considered (Burnham and Anderson 2002, p.68).
# #

# globx <- x <- mod

# chat <- summary(globx)$dispersion



# k <- length(coef(x)) + 1      # add 1 for the estimate of chat
# (deviance(x) / chat) + 2 * k

# x <- mod

# sum(x $ y * log(x$fitted.values) - x$fitted.values - lgamma(x$y + 1))



rm(list=ls())

# END TESTS
