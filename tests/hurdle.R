suppressPackageStartupMessages(library(pscl))
library(MuMIn)
data(bioChemists)

fmh1 <- hurdle(art ~ fem + mar*phd, data = bioChemists, dist = "negbin",
    zero.dist = "negbin")
(dd <- dredge(fmh1, m.min = 1))


#fm_pois <- glm(art ~ ., data = bioChemists, family = poisson)
#fm_qpois <- glm(art ~ ., data = bioChemists, family = quasipoisson)
#fm_nb <- glm.nb(art ~ ., data = bioChemists)
#dredge(fm_nb, m.max = 2)

## with simple inflation (no regressors for zero component)
fm_zip <- zeroinfl(art ~ . | 1, data = bioChemists)

(dd1 <- dredge(fm_zip, m.max = 2, eval = TRUE))
summary(ma <- model.avg(dd))
coef(ma)
confint(ma)

fm_zinb <- zeroinfl(art ~ . | phd + mar, data = bioChemists, dist = "negbin")
(dd2 <- dredge(fm_zinb, m.max = 1, eval = TRUE))
summary(model.avg(dd2))

fm_zinb <- zeroinfl(art ~ phd + mar, data = bioChemists, dist = "negbin")

(dd2 <- dredge(fm_zinb, m.max = 2, eval = TRUE, trace = T))

fmh1dot <- hurdle(art ~ ., data = bioChemists, dist = "negbin",
    zero.dist = "negbin")

dd2 <- dredge(fmh1dot, m.min = 1, m.max = 1, eval = T, trace = T)

mod <- get.models(dd2[1:4])

summary(model.avg(mod))
summary(model.avg(dd2[1:4]))
#summary(model.avg(dd2))
#coef(fmh1dot)
#coefTable(fmh1dot)
#summary(fmh1dot)
