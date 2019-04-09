if (MuMIn:::testStart("nnet", "MASS")) {

    # Trimmed-down model from example(birthwt)
    data(birthwt)

    bwt <- with(birthwt, data.frame(
            low = low,
            race = factor(race, labels = c("white", "black", "other")),
            ptd = factor(ptl > 0),
            smoke = (smoke > 0)
            ))

    options(contrasts = c("contr.treatment", "contr.poly"))
    bwt.mu <- multinom(low ~ ., data = bwt)
    dd <- dredge(bwt.mu, trace=T)

    summary(model.avg(dd[1:5]))
    gm <- get.models(dd, subset = 1:5)
    ma <- model.avg(gm)

    summary(ma)

    # predict(ma) // Cannot average factors!
}

