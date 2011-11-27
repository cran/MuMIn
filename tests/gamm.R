library(mgcv)
library(gamm4)
library(MuMIn)


RNGkind("Mersenne")
set.seed(0) # 16
dat <- gamSim(6, n=100, scale=5, dist="normal")

fmgs2 <- gamm(y ~s(x0)+ s(x3) + s(x2), family=gaussian, data=dat, random = list(fac=~1))

dd <- dredge(fmgs2)

#dd <- dredge(mgcv::gamm(y ~s(x0)+ s(x3) + s(x2), family=gaussian, data=dat, random = list(fac=~1)))
#dd <- dredge(gamm4::gamm4(y ~s(x0)+ s(x3) + s(x2), family=gaussian, data=dat,
	#random = ~(1|fac)), trace=T)
#mod <- gamm4::gamm4(y ~s(x0)+ s(x3) + s(x2), family=gaussian, data=dat, random = ~(1|fac))
#models <- get.models(dd, delta <= 4))
summary(model.avg(dd, delta <= 4))

#plot(dd, col=2, col2=8, col.bias=1)



fmg4s2 <- gamm(y ~s(x0)+ s(x3) + s(x2), family=gaussian, data=dat, random = ~ (1|fac))
dd4 <- dredge(fmg4s2)
summary(avg <- model.avg(dd4, delta <= 4))
