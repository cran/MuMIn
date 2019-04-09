if (MuMIn:::testStart("MASS")) {

    quine.nb1 <- glm.nb(Days ~ 0 + Sex / (Age + Eth * Lrn), data = quine)
    #quine.nb1 <- glm.nb(Days ~ Sex/(Age + Eth*Lrn), data = quine)

    ms <- dredge(quine.nb1)

    models <- get.models(ms, subset = TRUE)
    models <- get.models(ms, subset = NA)

    print(summary(model.avg(models)))

    #dredge(quine.nb1) # OK
    #dredge(quine.nb1x = NA) # OK
    #dredge(quine.nb1) # OK
    print(dredge(quine.nb1)) # OK
    #dredge(quine.nb1) # Right, should be the same as above
    ma <- model.avg(dredge(quine.nb1), subset = cumsum(weight) <= .9999)

    print(summary(ma))
    # Cannot predict with this 'averaging'
    #pred <- predict(ma, se=TRUE)

    #pred <- cbind(pred$fit, pred$fit - (2 * pred$se.fit), pred$fit + (2 * pred$se.fit))
    #matplot(pred, type="l")
    #matplot(family(quine.nb1)$linkinv(pred), type="l")
}