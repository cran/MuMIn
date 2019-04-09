if (MuMIn:::testStart("MASS")) {

    data(Cement, package = "MuMIn")

    nseq <- function(x, len = length(x)) 
        seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length = len)

    fm1 <- rlm(y ~X1 + X2 * X3 + X4, data = Cement)
    dd <- dredge(fm1, trace = TRUE)
    gm <- get.models(dd, subset = 1:10)
    ma <- model.avg(gm)
    stopifnot(all(predict(ma) == predict(ma, Cement)))
    predict(ma, lapply(Cement, nseq, len = 30), se.fit = TRUE)
    vcov(ma)
}