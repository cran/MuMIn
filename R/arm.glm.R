arm.glm <-
function(object, R = 250, weight.by = c("aic", "loglik"), trace = FALSE) {
	if(!inherits(object, c("glm", "lm")))
		stop("'object' must be a \"glm\" or \"lm\" object")
		
	maxtrials <- 10L
	weight.by <- switch(match.arg(weight.by), aic = 1L, loglik = 2L)
	allterms <- getAllTerms(object)
	ordtrm <- attr(allterms, "order")
	deps <- attr(allterms, "deps")[ordtrm, ordtrm]
	nterms <- length(allterms)
	mm <- model.matrix(object)
	n1 <- ceiling((nall <- nrow(mm))/2)
	n2 <- nall - n1 

	combinTerms <- lapply(seq.int(0L, 2^nterms - 1L), function(j)
		as.logical(intToBits(j)[1L:nterms]))
	combinTerms <- combinTerms[vapply(combinTerms,
		formula_margin_check, FALSE, deps)]
	## combinPredictors
	combin <- vapply(combinTerms, function(x, assign)
		assign %in% c(0, which(x)), logical(ncol(mm)),
		assign = attr(mm, "assign"))
	combin <- combin[, i <- colSums(combin) < min(n1, n2) - 1L, drop = FALSE]
	combinTerms <- combinTerms[i]
	
	nModels <- ncol(combin)
	fam <- family(object)
	y <- object$y
	off <- object$offset
	if(is.null(y)) y <- get.response(object)

	prwt <- weights(object, "prior")
	if(is.null(prwt)) prwt <- rep(1, nall) 
	
	if(fam$family == 'binomial') {
		y <- prwt * cbind(y, 1 - y, deparse.level = 0L)
		yvectorize <- function(y) y[, 1L] / rowSums(y)
		prwt <- rep(1, NROW(y))
		jresp <- c(1L, 2L)
	} else {
		yvectorize <- function(y) y[, 1L] ## XXX === drop  ?
		jresp <- 1L
	}
	Z <- cbind(y, mm, prwt, deparse.level = 0L)
	rownames(Z) <- NULL
	jprwts <- ncol(Z)
	jterms <- (1L:(ncol(Z) - 1L))[-jresp]
	
	traceinfo <- if(trace) { function(...) cat(..., sep = "") } else function(...) {}
	
	## begin iteration
	wts <- matrix(NA_real_, ncol = nModels, nrow = R)
	for(iter in 1L:R) {
		traceinfo("iteration=", iter, "\n")
		for (u in seq.int(maxtrials)) {
			traceinfo("  trial=", u, "\n    ")
			i <- sample.int(nall, n1)
			y1 <- Z[i, jresp, drop = FALSE]
			y2 <- Z[-i, jresp, drop = FALSE]
			vy2 <- yvectorize(y2)
			x1 <- Z[i, jterms, drop = FALSE]
			x2 <- Z[-i, jterms, drop = FALSE]
			prwts1 <- Z[i, jprwts]
			prwts2 <- Z[-i, jprwts]
			ic <- numeric(nModels)
			for (k in seq.int(nModels)) {
				traceinfo("k=", k, " ")
				fit1 <- glm.fit(x1[, combin[, k], drop = FALSE], y1, family = fam, weights = prwts1,
								offset = off)
				if (problem <- (anyNA(fit1$coefficients) || !fit1$converged)) {
					traceinfo("<!>")
					break
				}
				ic[k] <- aicloglik_glm_fit(fit1, vy2, x2[, combin[, k], drop = FALSE], 
					prwts2, off)[weight.by]
			}
			traceinfo("\n")
			if (!problem) 
				break
		}
		d <- exp(-ic / 2)
		wts[iter, ] <- d/sum(d)
	}
	
	
	wts[!is.finite(wts)] <- NA_real_
	wts <- wts[round(rowSums(wts, na.rm = TRUE)) == 1, ]
	wts <- wts/rowSums(wts, na.rm = TRUE)
	wtsmean <- colMeans(wts)
	
	cfnames <- colnames(Z[, jterms])
	x <- Z[, jterms, drop = FALSE]
	y <- Z[, jresp, drop = FALSE]
	msTable <- matrix(NA_real_, nrow = nModels, ncol = 6L, dimnames = list(1L:nModels, 
		c("df", "logLik", "AIC", "delta", "weight", "ARM weight")))
	## coefArray(model, c(coef, se, df), coefficients)
	coefArray <- array(NA_real_, dim = c(nModels, 3L, length(cfnames)), dimnames = list(1L:nModels, 
		c("Estimate", "Std. Error", "df"), cfnames))
	for (k in seq.int(nModels)) {
		fit1 <- glm.fit(x[, combin[, k], drop = FALSE], y, family = fam, offset = off)
		coefArray[k, , combin[, k]] <- rbind(t(summary.glm(fit1)$coefficients[, 1L:2L, 
			drop = FALSE]), fit1$df.residual)
		msTable[k, c(3L, 2L, 1L)] <- aicloglik_glm_fit(fit1, fit1$y, x[, combin[, k], 
			drop = FALSE], fit1$prior.weights, off)
	}
	msTable[, 4L] <- msTable[, 3L] - min(msTable[, 3L])
	msTable[, 5L] <- Weights(msTable[, 3L])
	msTable[, 6L] <- wtsmean
	
	cfmat <- coefArray[, 1L, ]
	cfmat[is.na(cfmat)]<- 0
	coefMat <- array(dim = c(2L, ncol(cfmat)),
		dimnames = list(c("full", "subset"), colnames(cfmat)))
	coefMat[1L, ] <- drop(wtsmean %*% cfmat)
	
	#debug <- list(wtsmean = wtsmean, cfmat = cfmat, coefArray = coefArray)
	
	ass <- attr(mm, "assign")
	bp <- !is.na(coefArray[, 1L, ass != 0L & !duplicated(ass)])
	tenm <- allterms[ordtrm]
	allmodelnames <- .modelNames(allTerms = apply(bp, 1L, function(z) tenm[z]),
						uqTerms = tenm)
	rownames(msTable) <- c(allmodelnames)
	
	ordmod <- order(msTable[,4L], decreasing = FALSE)
		
	rval <- list(
		msTable = structure(as.data.frame(msTable[ordmod, ]),
			term.codes = attr(allmodelnames, "variables")),
		coefficients = coefMat,
		coefArray = coefArray[ordmod, , ],
		sw = {
			structure(wtsmean %*% bp,
				n.models = structure(colSums(bp), names = tenm),
				names = tenm, class = "sw")
		 },
		formula = object$formula,
		call = match.call() 
	    #, debug = debug
	)
		
	attr(rval, "rank") <- .getRank(AIC)  ## TODO
	attr(rval, "nobs") <- nrow(x)
	attr(rval, "beta") <- "none"
	attr(rval, "revised.var") <- TRUE
	attr(rval, "arm") <- TRUE
	attr(rval, "model.weights") <- "ARM"

	class(rval) <- "averaging"
	rval 
}


armWeights <-
function(object, ..., data, weight.by = c("aic", "loglik"), R = 1000) {
	weight.by <- switch(match.arg(weight.by), aic = 1L, loglik = 2L)
    
    models <- getModelArgs()
    m <- length(models)
    if(m < 2) stop("need more than one model")
		
	p <- 0.5
		
	if(!all(vapply(models, inherits, TRUE, "lm")))
	   stop("'models' must inherit from \"lm\" class")
  
    R <- as.integer(R[1L])
    if(R <= 1) stop("'R' must be positive")

    n <- nrow(data)
    nt <- ceiling(n * p)
	
	maxtrials <- 3L
    wmat <- array(dim = c(R, m))
	r <- counter <- 1L
    counterLimit <- R * maxtrials # 	
    mode(R) <- mode(counterLimit) <- "integer"
    while(counter < counterLimit && r <= R) {
        counter <- counter + 1L
        k <- sample.int(n, size = nt)
		data.test <- data[-k, , drop = FALSE]
		data.train <- data[k, , drop = FALSE]
		y.test <- get.response(models[[1L]], data.test)
        for(j in seq.int(m)) {
			fit <- models[[j]]
            tf <- terms(fit)
            fam <- family(fit)
			off <- fit$offset
			if(is.null(off)) off <- rep(0, n)
            wts <- fit$weights
			if(is.null(wts)) wts <- rep(1, n)
			fit1 <- do_glm_fit(tf, data.train, family = fam, weights = wts[k], offset = off[k])
			if(!fit1$converged) break
			wmat[r, j] <- aicloglik_glm_fit(fit1, y.test, model.matrix(tf, data.test), wts[-k], off[-k])[weight.by]
        }
		if(!anyNA(wmat[r, ]))
			r <- r + 1L
    }
	
	wmat <- exp(-wmat / 2)
    wts <- colMeans(wmat)
    wts <- wts / sum(wts)
	structure(wts, wt.type = "ARM", names = names(models), class = c("model.weights", "numeric"))
}

