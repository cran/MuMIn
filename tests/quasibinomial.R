library("MuMIn")
options(na.action = na.fail)

budworm <- data.frame(ldose = rep(0:5, 2), numdead = c(1, 4, 9, 13, 18, 20, 0,
	2, 6, 10, 12, 16), sex = factor(rep(c("M", "F"), c(6, 6))))
budworm$SF = cbind(numdead = budworm$numdead, numalive = 20 - budworm$numdead)

qbinomial <- function(...) {
	res <- quasibinomial(...)
	res$aic <- binomial(...)$aic
	res
}

budworm.qqlg <- glm(SF ~ sex*ldose + sex*I(ldose^2), data = budworm, family = qbinomial)
budworm.lg <- glm(SF ~ sex*ldose + sex*I(ldose^2), data = budworm, family = binomial)

r.squaredLR(budworm.lg)
r.squaredGLMM(budworm.lg)

dd <- dredge(budworm.lg, rank = "QAIC", chat = summary(budworm.lg)$dispersion)
#dd <- dredge(budworm.lg) # should be the same
mod <- get.models(dd, subset = NA)
mod

# Note: this works:
# ma <- model.avg(mod)
# but this will not: ('rank' attribute passed from 'dredge' is lost)
# ma <- model.avg(mod)
# so, need to supply them

model.avg(mod[1:5], rank = "QAICc", rank.args = list(chat = 0.403111))


