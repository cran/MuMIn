predict.averaging.multinom <-
function(object, ...,
        predict = function(model) as.matrix(predict(model, type = "prob", ...))
        ) {
    models <- get.models(object, TRUE)
    pp <- lapply(models, predict, ...)
    pa <- array(dim = c(length(pp), dim(pp[[1L]])))
    for(i in seq_along(models))
        pa[i, , ] <- pp[[i]]
    # model-averaged predicted class probabilities:    
    pavg <- apply(pa, c(2L, 3L), weighted.mean, weights = Weights(object))
    models[[1L]]$lev[max.col(pavg)]
}