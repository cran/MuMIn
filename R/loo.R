
loo <-
function(object, type = c("loglik", "rmse"), ...)
	UseMethod("loo")

loo.default <-
function(object, type = c("loglik", "rmse"), ...)
	.NotYetImplemented()
	# for other types of models use manipulated call - SLOW

# lm, glm, mgcv::gam
loo.lm <- 
function(object, type = c("loglik", "rmse"), start, etastart, mustart,
		 control, intercept, ...) {
	
	if(!inherits(object, "lm"))
		stop("'object' must be a \"glm\" or \"lm\" object")
	
	## TODO: pass other arguments: 
	if(!missing(start) || !missing(etastart)  || !missing(mustart)  ||
	   !missing(control) || !missing(intercept)) {
		warning("arguments 'start', 'etastart', 'mustart', 'control', 'intercept' are ignored")
	}
	## binary response: object$y is always a vector and weights(object) give size

	if(type != "rmse_t") type <- match.arg(type) # hidden option
	tt <- terms(object)
	beta <- coef(object)
	fam <- family(object)
	X <- model.matrix(object)
	y0 <- get.response(object) # model.frame(object)[, asChar(attr(tt, "variables")[-1L][[attr(tt, "response")]])]
	nobs <- NROW(y0)
	wt <- weights(object, "prior")
	if(is.null(wt)) wt <- rep(1, nobs)
	
	if(NCOL(y0) == 2L) { # binomial
		# NOTE: don't need to keep 'y' as 2-column matrix, but if also weights
		# 		are given, 'glm.fit' warns about "non-integer # of successes"
		n <- rowSums(y0)
		y <- y0[, 1L] / n
		wt0 <- wt / n
		# FUNFACT:  'x[1L] == x' more efficient than x[1L] == x[-1L]
		# If weights are equal, use y as proportion, otherwise a matrix
		if(all(wt == n)) { #if(all(wt0[1L] == wt0)) {
			y0 <- y
			wt0 <- wt
		}
	} else  {
		n <- rep(1, nobs)
		y <- y0
		wt0 <- wt
	}
	y0 <- as.matrix(y0)

	offset <- object$offset
	if(is.null(offset)) offset <- numeric(nobs)
	
	func <- switch(type,
	loglik = {
		dev <- function(y, mu, wt, fam)  sum(fam$dev.resids(y, mu, wt))
		llik <- function(y, X, beta, fam, n, wt = 1, off = NULL) {
			# wt : fit$prior.weights
			wt <- rep(wt, length.out = NROW(y))
			mu <- predict_glm_fit(beta, X, off, fam)[, 1L]
			ep <- if (fam$family %in% c("gaussian", "Gamma", "inverse.gaussian")) 1 else 0
			(fam$aic(y, n, mu, wt, dev(y, mu, wt, fam)) / 2) - ep # +LL
		}	
		function(fit, z) {
			llik(z[["y"]], z[["X"]], fit$coefficients, fit$family, z[["n"]],
				z[["wt"]], z[["offset"]])
			}
		},
	rmse = function(fit, z) {  # RMSE on response scale
		beta <- fit$coefficients
		py <- predict_glm_fit(beta, z[["X"]], z[["offset"]], fit$family)[, 1L]
		z[["y"]] - py # inefficient to '^2' here, do it later
	},
	rmse_t = function(fit, z) { # alternatively: MSE on transformed data
		beta <- fit$coefficients
		py <- predict_glm_fit(beta, z[["X"]], z[["offset"]])[, 1L]
	    # prediction on link scale
		z[["eta"]] - py # inefficient to '^2' here, do it later
		# XXX: problem with RMSE with binomial and y = 0 or 1 (= Inf on link scale) 
	})
	
	.DebugPrint(type)
	.Debug(if(type == "loglik") {
		.DebugPrint(y0)
		message("running test 1...")
		# XXX: DEBUG test // lm has no $family
		testLL1 <- llik(y, X, object$coefficients, family(object), n, wt, offset)
		
		.DebugPrint(testLL1)
		.DebugPrint(logLik(object))
		.DebugPrint(testLL1 - logLik(object))
		.DebugPrint(rbind(n, wt, offset))

		stopifnot(all.equal(-testLL1, c(logLik(object)), tolerance = 1e-5))
		message("OK")
		message("running test 2...")
		testFm <- glm.fit(y = y0, x = X, family = fam, offset = offset, weights = wt0)
		.DebugPrint(rbind(testFm$coefficients, object$coefficients))
		stopifnot(all.equal(testFm$coefficients, object$coefficients))
		message("OK")
		message("running test 3...")
		testLL2 <- llik(y, X, testFm$coefficients, testFm$family, n, wt, offset)
		.DebugPrint(c(testLL2, logLik(object)))
		stopifnot(all.equal(-testLL2, c(logLik(object)), tolerance = 1e-5))
		message("OK")
		.DebugPrint(testLL2)
		.DebugPrint(logLik(object))
	})
	
	fdat <- if(type == "loglik") {
		data.frame(X = 0, y = y, n = n, wt = wt, offset = offset)
	} else {
		data.frame(X = 0, y = y, eta = fam$linkfun(y), offset = offset)
	}
	fdat$X <- X
	
	rval <- numeric(nobs)
	for (i in seq.int(nobs)) {
		fm1 <- glm.fit(y = y0[-i, , drop = FALSE], x = X[-i, , drop = FALSE],
				family = fam, offset = offset[-i], weights = wt[-i])
		rval[i] <- func(fm1, fdat[i, ])
	}
	
	if(type == "loglik") mean(rval) else sqrt(mean(rval^2))
}
