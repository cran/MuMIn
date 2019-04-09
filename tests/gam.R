if (MuMIn:::testStart("mgcv")) {

    RNGkind("Mersenne")
    set.seed(0) ## simulate some data...
    dat <- gamSim(1, n = 400, dist = "binary", scale = 2)
    #gam1 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat)

    ops <- options(warn = -1)

    gam1 <- gam(y ~ s(x0) + s(x1) + s(x2) +  s(x3) + (x1+x2+x3)^2,
        data = dat, method = "GCV.Cp", family = binomial)

    dd <- dredge(gam1, subset = !`s(x0)` & (!`s(x1)` | !x1) & 
        (!`s(x2)` | !x2) & (!`s(x3)` | !x3), fixed = "x1")
        
    gm <- get.models(dd, cumsum(weight) <= .95)
    ma <- model.avg(gm)

    print(summary(ma))

    print(predict(ma, dat[1:10, ], se.fit = TRUE, type = "link"))
    print(predict(ma, dat[1:10, ], se.fit = TRUE, type = "response"))
    print(predict(ma, dat[1:10, ], se.fit = TRUE, type = "link", 
        backtransform = TRUE))

    options(ops)
}
