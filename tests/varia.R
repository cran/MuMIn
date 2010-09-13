# Test varia

require(MuMIn)

print(packageDescription("MuMIn", fields = "Version"))


# TEST binary response --------------------------------------------------------------
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

# END TESTS