

# TODO: common framework for glm.fit loops
# TODO: pass data in formula environment

BGWeights <-
function(object, ..., data, force.update = FALSE) {
     
    models <- getModelArgs()
    M <- length(models)
    if(M < 2) stop("need more than one model")
    checkIsModelDataIdentical(models)

    no <- nrow(data)
    k <- sample.int(no, floor(no / 2))
    dat_train <- data[k, ]
    dat_test <- data[-k, ]
    # TODO: allow user to specify offset and weights  
    offset <- rep(0, no) 
    weights <- rep(1, no)
    weights_train <- weights[k]
    weights_test <- weights[-k]
    offset_train <- offset[k]
    offset_test <-  offset[-k]
    
    if(!force.update && all(vapply(models, inherits, FALSE, "glm"))) { # XXX: what about lm
        # XXX: here 'offset_train' DOES include offset specified in 'formula'
        
        py_test <- array(dim = c(nrow(dat_test), M))
        for(i in seq.int(length(models))) {
            fit <- models[[i]]
            tf <- terms(fit)
            fit_train <- update_glm_fit(fit, dat_train, weights_train, offset_train)
            py_test[, i] <- predict_glm_fit(fit_train$coefficients, model.matrix(tf, dat_test),
                offset = offset_test, family = family(fit))[, 1L]
        }
       
    } else { # for non-glm models, use update (2x slower)
        
        if(any(!vapply(models, function(x) is.null(get_call(x)$offset), FALSE)))
            stop("use of 'offset' argument in model calls is not allowed. ",
                 "Specify 'offset' term in the formula instead")
            
        # NOTE: better to update in `parent.frame` because of possible masking
        #       of the variables in model's call by variables defined in this
        #       function. For example:
        ##      fm <- glm(y ~ X1, data = Cement)
        ##      (function(x) {
        ##          Cement <- NA
        ##          update(x) # Error
        ##      })(fm)
        ##      Even better is to provide all variables in environment(formula)
        ##      and rename the formula variables
        pf <- parent.frame()
        
        z_data <- tmpvarname(pf)
        z_weights <- tmpvarname(pf)
        assign(z_data, dat_train, pf)
        assign(z_weights, weights_train, pf)
        
        on.exit(rm(list = c(z_data, z_weights), envir = pf, inherits = FALSE))
        nz_data <- as.name(z_data)
        nz_weights <- as.name(z_weights)
        
        # XXX: here 'offset_train' DOES NOT include offset specified in 'formula'
        py_test <- array(dim = c(nrow(dat_test), M))
        for(i in seq.int(M)) {
            cl <- get_call(models[[i]])
            cl$data <- nz_data
            cl$weights <- nz_weights
            #print(cl)
            py_test[, i] <- predict(eval(cl, pf), newdata = dat_test,
                 type = "response")
                # XXX: weird behaviour of predict when offset= is given.
                #      prediction always has length of the offset. a bug in predict.glm?
                #offset = `*tmp_dat*`[[3L]]),
        }
        # NOTE: No speed gain with *apply - more memory needed to store train_fits
        ## train_fits <- lapply(models, function(x) update(x, data = `*tmp_dat`[[1L]], weights = `*tmp_dat`[[2L]], offset = `*tmp_dat`[[3L]]))
        ## py_test <- vapply(train_fits, predict, numeric(nrow(dat_test)), newdata = dat_test, type = "response")

    }
    
    y_test <- get.response(models[[1L]], data = dat_test)
    if(is.matrix(y_test)) y_test <- y_test[, 1L] / rowSums(y_test) # binomial
    
    Sigma <- cov(y_test - py_test)
    # XXX: I want do avoid dependency on MASS
    #ginv <- if(use.MASS) getFrom("MASS", "ginv") else solve
    
    ones <- rep(1, M)
    
    fn1 <- function(ones, Sigma, ginv) ginv(t(ones) %*% ginv(Sigma) %*% ones) %*% ones %*% ginv(Sigma)

    rval <- tryCatch(fn1(ones, Sigma, solve), error = function(e) {
        if(length(find.package("MASS", quiet = TRUE)) == 1L)
            fn1(ones, Sigma, getFrom("MASS", "ginv")) else
            stop(e)
        })[1L, ]
    
    structure(rval, wt.type = "Bates-Granger", names = names(models),
        class = c("model.weights", class(rval)))
}

