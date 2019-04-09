library("MuMIn")
options(na.action = "na.fail")
data(Orthodont, package = "nlme")
    
fm1 <- lm(distance ~ Sex * age + age * Sex, data = Orthodont)

dispersion <- function(object) {
    wts <- weights(object)
    if (is.null(wts)) 
        wts <- 1
    sum((wts * resid(object, type = "working")^2)[wts > 0])/df.residual(object)
}

dd <- dredge(fm1, extra = alist(dispersion))
gm <- get.models(dd, subset = 1:4)
ma <- model.avg(gm, revised = F)

vcov(ma)
summary(ma)
confint(ma)

predict(ma)
predict(ma, se.fit = TRUE)
predict(ma, data.frame(Sex = "Male", age = 8:12))

rm(list = ls())

data(Cement, package = "MuMIn")

nseq <- function(x, len = length(x)) seq(min(x, na.rm = TRUE), 
    max(x, na.rm = TRUE), length = len)

fm1 <- glm(y ~ (X1 + X2 + X3)^2, data = Cement)
dd <- dredge(fm1)

gm <- get.models(dd, subset = 1L:10L)

summary(ma <- model.avg(gm))

vcov(ma)

summary(ma1 <- model.avg(dd[1L:10L]))
summary(ma2 <- model.avg(model.sel(dd[1L:10L], rank = "AICc")))
all.equal(ma$avg.model, ma1$avg.mode)

predict(ma) == predict(ma, Cement)
predict(ma, se.fit = TRUE)
predict(ma, lapply(Cement, nseq))
