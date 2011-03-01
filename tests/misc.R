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
budworm$SF = cbind(numdead = budworm$numdead, numalive = 20 - budworm$numdead)


(function(dat) (function(dat2) {
	mod <- glm(SF ~ sex*ldose, data = dat2, family = "quasibinomial")
	dd <- dredge(mod, rank = "QAIC", chat = summary(budworm.lg)$dispersion)
	gm <- get.models(dd, family="binomial")
	#print(sys.frames())
	model.avg(gm)
})(dat))(budworm)

rm(list=ls())

# END TESTS
