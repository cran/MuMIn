


jackknifeWeights <- 
function(object, ...,
		 data,
		 type = c("loglik", "rmse"),
		 family = NULL, weights = NULL,
		 optim.method = "BFGS",
		 maxit = 1000,
		 # seed = NULL,
		 optim.args = list(),
		 start = NULL,
		 force.update = FALSE,
		 py.matrix = FALSE
		) {
	
	models <- getModelArgs()
	M <- length(models)
	if(M < 2L) stop("need more than one model")
	checkIsModelDataIdentical(models)
	
	type <- match.arg(type)[1L]
	useWeightsArg <- !missing(weights) && !is.null(weights)
	hasPymat <- is.matrix(py.matrix)
	
	if(!is.null(start)) {
		if(length(start) != M)
			stop("length of 'start' should equal the number of models [", M, "]")
		start <- log(start)
		start <- (start - start[1L])[1L]
	} else {
		start <- runif(M - 1L)
	}

	if(type == "loglik") {
		if(is.null(family)) {
			famstr <- tryCatch(vapply(models, function(m) {
					f <- family(m)
					c(f[["link"]], f[["family"]])
				}, character(2L)), error = function(e) NULL)
				if (is.null(famstr))
					stop("cannot get 'family' function for some models")
				if (!all(famstr[1L, ] == famstr[1L, 1L]) && all(famstr[2L, ] == famstr[2L, 1L]))
					stop("models use different 'family' and 'family' argument is not specified")
			family <- stats::family(models[[1L]])
		} else {
			# code from 'stats::glm'
			if (is.character(family)) 
				family <- get(family, mode = "function", envir = parent.frame())
			if (is.function(family)) 
				family <- family()
			if (is.null(family$family)) {
				stop("'family' not recognized")
			}
		}
	}
	
	# generate matrix of predictions
	if(hasPymat) {
		pymat <- py.matrix
		Y <- as.matrix(get.response(models[[1L]], if(missing(data)) NULL else data))
	} else if(!force.update && all(vapply(models, inherits, FALSE, "glm"))) {
		
		if (isTRUE(getOption("debug.MuMIn"))) message("using glm.fit")
		
		mf <- mergeMF(models, check = FALSE)
		tf <- terms(mf)
		if(!missing(data)) mf <- model.frame(tf, data = data)
	
		X <- model.matrix(tf, mf)
		Y <- as.matrix(get.response(mf))
		no <- NROW(Y)

		colnames(X) <- fixCoefNames(colnames(X))
		xil <- lapply(models, function(fit, allcn) {
			match(fixCoefNames(names(coef(fit))), allcn)
			}, colnames(X))

		fitting_wts <- numeric(no)
		fitting_wts[] <- if(useWeightsArg) weights else 1

		pymat <- array(dim = c(no, M))
		xseq <- seq.int(no)
		for(j in seq.int(M)) {
			fit <- models[[j]]
			tf <- terms(fit)
			fam <- family(fit)
			off <- fit$offset
			
			if(!useWeightsArg) {
				fitting_wts <- fit$prior.weights
				if(fam$family == "binomial" && NCOL(Y) == 2L)
					fitting_wts <- fitting_wts / rowSums(Y)
			}
			if(is.null(off)) off <- rep(0, NROW(Y))
				
			.Debug({# DEBUG
				message("testing glm.fit #",  j)
				cf1 <- glm.fit(X[, xil[[j]], drop = FALSE], Y, family = family(fit),
					weights = fitting_wts, offset = off)$coefficients
				cf2 <- models[[j]]$coefficients
				names(cf2) <- fixCoefNames(names(cf2))
				stopifnot(all.equal(cf1[names(cf2)], cf2))
			}) # DEBUG
		
			for(i in xseq) {
				coef1 <- glm.fit(X[-i, xil[[j]], drop = FALSE], Y[-i, , drop = FALSE],
					family = family(fit),
					weights = fitting_wts[-i], offset = off[-i])$coefficients
				pymat[i, j] <- predict_glm_fit(coef1, X[i, xil[[j]], drop = FALSE],
					offset = off[i], family = fam)[, 1L]
			}
		}
	} else {
		if(isTRUE(getOption("debug.MuMIn")))
			message("using update")
		
		if(any(!vapply(models, function(x) is.null(get_call(x)$offset), FALSE)))
			stop("use of 'offset' argument in model calls is not allowed. ",
				 "Specify 'offset' term in the model formula instead")
			
		Y <- as.matrix(get.response(models[[1L]], data))
		no <- NROW(Y)
		# TODO: check for prior weights consistency between models
		fitting_wts <- numeric(no)
		fitting_wts[] <- if(useWeightsArg) weights else 1
		if(!useWeightsArg) {
			fitting_wts <- weights(models[[1L]], "prior")
			if(is.null(fitting_wts))
				fitting_wts <- rep(1, no) else
				if(NCOL(Y) == 2L)
					fitting_wts <- fitting_wts / rowSums(Y)
		}

		pf <- parent.frame()
		
		# need some variable names unused in parent.frame
		z_data <- tmpvarname(pf)
		z_subset <- tmpvarname(pf)

		cll <- eval(substitute(lapply(models, function(x) {
			update(x, evaluate = FALSE,
				data = DATA[[1L]][INDEX, ],
				weights = DATA[[2L]][INDEX]
				)}), list(DATA = as.name(z_data), INDEX = as.name(z_subset))))
		
		assign(z_data, list(data, fitting_wts), pf)
		on.exit(rm(list = c(z_data, z_subset), envir = pf, inherits = FALSE))
		
		pymat <- array(dim = c(no, M))
		for(i in 1L:no) {
			assign(z_subset, -i, pf)
			pymat[i, ] <- vapply(cll,
				function(x, newdata, envir) predict(eval(x, envir),
					newdata = newdata, type = "response"),
				newdata = data[i, , drop = FALSE],
				envir = pf, FUN.VALUE = numeric(1L))
		}
	} ## produces pymat
	
	if(isTRUE(py.matrix)) {
		dimnames(pymat) <- list(NULL, names(models))
		return(pymat)
	}
	
	if("control" %in% names(optim.args)) {
		optim.control <- optim.args$control
		if(is.numeric(optim.args$control$maxit) && missing(maxit))
			maxit <- optim.args$control$maxit
		optim.args$control <- NULL
	} else optim.control <- list()
	optim.control$maxit <- maxit
	
	if(NCOL(Y) == 2) {
		ns  <- rowSums(Y)
		y <- Y[, 1L] / ns
	} else {
		no <- NROW(Y)
		y <- c(Y)
		ns <- rep(1, no)
	}
	
	if("loglik" == type) {
	#if("loglik" %in% type) {
		llik <-
		function(y, mu, fam, n, wt) {
			wt <- rep(wt, length.out = NROW(y)) # wt : fit$prior.weights
			ep <- if (fam$family %in% c("gaussian", "Gamma", "inverse.gaussian")) 1 else 0
			(fam$aic(y, n, mu, wt, sum(fam$dev.resids(y, mu, wt))) / 2) - ep # +LL
		}

		if(useWeightsArg) {
			prior_wts <- if(NCOL(Y) == 2) {
				fitting_wts * ns
			} else fitting_wts
		} else
			prior_wts <- weights(models[[1L]], "prior")
		
		optfun_loglik <- 
		function(w1, y, pymat, fam, n, wt) {
			w <- c(1, exp(w1))
			w <- w / sum(w)
			llik(y, pymat %*% w, fam, n, wt)
		}
	
	
		optres <- eval(as.call(c(alist(optim, optfun_loglik,
				y = y, pymat = pymat, fam = family, n = ns, wt = prior_wts,
				par = start,
				method = optim.method,
				control = optim.control),
			optim.args)))
		
	} else { #if("rmse" == type) {
	#if("rmse" %in% type) {
		# compute RMSE for a value of w, given J: 
		optfun_rmse <- function(ww, pymat, y) {
			w <- c(1, exp(ww))
			w <- w / sum(w)
			py <- pymat %*% w
			sqrt(mean((py - y)^2))
		}
		optres <- eval(as.call(c(alist(optim, optfun_rmse,
				pymat = pymat, y = y,
				par = runif(M - 1L),
				method = optim.method, control = optim.control),
			optim.args)))
	
	}
	
	if (optres$convergence != 0)
		stop("not converged. 'optim' gave error code [", optres$convergence, "]")

	wts <- exp(optres$par)	
	wts <- c(1, wts)
	wts <- wts / sum(wts)
	
	structure(wts, wt.type = paste0("jackknife [", type, "]"),
			  names = names(models),
			  class = c("model.weights", class(wts)))
} # end jackknifeWeights



#get_response_size_weigths <-
#function(y, w = rep(1, NROW(y))) {
#	if(NCOL(y) == 2L) {
#		n <- rowSums(y)
#		y <- y[, 1L] / n
#		w <- w / n
#	} else {
#		n <- rep(1, NROW(y))
#	}
#	list(y = y, n = n, weigths = w)
#}
