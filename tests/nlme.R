if (MuMIn:::testStart("nlme")) {
    fm1Dial.gls <- gls(rate ~ (pressure + I(pressure^2) + I(pressure^3)) * QB, Dialyzer, 
        method = "ML") 

    varying <- list(correlation = alist(AR1_0.771 = corAR1(0.771, form = ~1 | Subject), 
        AR1 = corAR1(), NULL), weights = alist(vp.press = varPower(form = ~pressure), 
        NULL)) 

    dd <- dredge(fm1Dial.gls, m.lim = c(1, 2), fixed = ~pressure, varying = varying)

    models <- get.models(dd, subset = 1:4)

    predict(fm1Dial.gls, se.fit = TRUE, newdata = Dialyzer[1:5, ])

    subset(dd, correlation == "AR1_0.771", recalc.delta = TRUE)

    ma <- model.avg(models, revised = TRUE)
    ms <- model.sel(models)
    print(ms, abbr = FALSE)
    print(ms, abbr = TRUE)
    summary(ma)
    predict(ma)[1:10] 

    # testing predict replacement:
    fm1 <- lme(rate ~ (pressure + I(pressure^2) + I(pressure^3)) * QB, ~1 | Subject, 
        data = Dialyzer)
    predict(fm1, newdata = Dialyzer[1:5, ], level = 0, se.fit = TRUE)
}     


