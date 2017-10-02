stackingWeights <-
function(object, ..., data, R, p = 0.5
         #, seed = NULL
         ) {
    
    models <- getModelArgs()
    M <- length(models)
    # TODO: check for lm class  / add support for other models
    if(M < 2L) stop("need more than one model")

    R <- as.integer(R[1L])
    if(R <= 1L) stop("'R' must be positive")
    if(p <= 0 || p >= 1) stop("'p' must be in range <0,1>")

    n <- nrow(data)
    nt <- round(n * p)

    wmat <- array(dim = c(R, M))
    r <- counter <- 1L
    counterLimit <- R * 2L
    mode(R) <- mode(counterLimit) <- "integer"
    while(counter < counterLimit && r <= R) {
        counter <- counter + 1L
        k <- sample.int(n, size = nt)
        pymat <- array(dim = c(n - nt, M))
        for(j in 1L:M) {
            fit <- models[[j]]
            tf <- terms(fit)
            fam <- family(fit)
            off <- fit$offset
            wts <- fit$weights
            coef1 <- do_glm_fit(tf, data[k, , drop = FALSE], fam, wts[k],
                off[k])$coefficients
            pymat[, j] <- predict_glm_fit(coef1,
                model.matrix(tf, data[-k, , drop = FALSE]), offset = off[-k],
                    family = fam)[, 1L]
        }
        y.test <- get.response(fit, data[-k, , drop = FALSE])
        sw1 <- tryCatch(.stacking(pymat, y.test), error = function(...) NULL)
        if(!is.null(sw1)) {
            wmat[r, ] <- sw1
            r <- r + 1L
        }
    }
    wts <- rbind(colMeans(wmat), apply(wmat, 2L, median),
                 deparse.level = 0)
    dimnames(wts) <- list(c("mean", "median"), names(models))
    
	structure(wts / rowSums(wts), wt.type = "stacking",
        class = c("model.weights", class(wts)))
}


.stacking <-
function(predicted, observed) {
   
    if(!is.matrix(predicted))
        stop("\"predicted\" must be a matrix")
    if(nrow(predicted) != length(observed))
        stop("number of rows in \"predicted\" is not equal to length of \"observed\"")

    if (NCOL(predicted) >= length(observed))
        stop("more models than test points. ",
             "Increase the test set or reduce the number of models")
        # TODO: make the error message more specific.

    # now do an internal splitting into "folds" data sets:
    weightsopt <-
    function(ww) {
        # function to compute RMSE on test data
        w <- c(1, exp(ww))
        w <- w / sum(w) 
        ## w all in (0,1) SIMON; set weight1 always to 1, other weights are
        ## scaled accordingly (this leads to a tiny dependence of optimal
        ## weights on whether model1 is any good or utter rubbish; see by moving
        ## the 1 to the end instead -> 3rd digit changes)
        pred <- as.vector(predicted %*% w)
        return(sqrt(mean((pred - observed)^2)))
    }

    ops <- optim(par = runif(NCOL(predicted) - 1L), weightsopt, method = "BFGS")
    if (ops$convergence != 0) stop("optimization not converged")
    round(c(1, exp(ops$par)) / sum(c(1, exp(ops$par))), 4L)
}


