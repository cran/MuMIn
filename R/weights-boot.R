


bootWeights <- 
function(object, ...,
		 R,
		 rank = c("AICc", "AIC", "BIC")
		) {
	
	models <- getModelArgs()
	M <- length(models)

	if(M < 2) stop("need more than one model")

	checkIsModelDataIdentical(models)
	
	for(fm in models) {
	  if(anyNA(match("x", names(fm)))) {
		warning("for efficiency of the bootstrap procedure, 'glm' should be called with 'x = TRUE'")
		break
	  }
	}
	
	rank <- match.arg(rank)
	
	ic <- switch(rank, AICc = AICc, AIC = AIC, BIC = BIC)
	
	mseq <- seq.int(M)
	
	if(all(vapply(models, inherits, FALSE, "glm"))) { ## !force.update &&
	  	best <- integer(R)
		ics <- numeric(M)
		no <- nobs(models[[1L]]) # assuming nobs is the same across models
		r <- 1L
		while(r <= R) {
			g <- sample.int(no, replace = TRUE)

			for(j in mseq) {
				fit <- models[[j]]
				fit1 <- glm.fit(model.matrix(fit)[g, , drop = FALSE], fit$y[g], fit$prior.weights[g],
				  offset = fit$offset[g], family = fit$family)
				ics[j] <- ic(loglik_glm_fit(fit1))
			}
			best[r] <- which.min(ics)
			r <- r + 1L
		}
	} else stop("all model objects must be of class \"glm\"")
  
	wts <- tabulate(best, M)
	wts <- wts  / sum(wts)
	structure(wts, wt.type = "bootstrap", names = names(models),
			  class = c("model.weights", class(wts)))
}




bootWeights2 <- 
function(object, ...,
		 R,
		 rank = c("AICc", "AIC", "BIC")
		) {
	
	models <- getModelArgs()
	M <- length(models)

	if(M < 2) stop("need more than one model")

	checkIsModelDataIdentical(models)
	
	for(fm in models) {
	  if(anyNA(match("x", names(fm)))) {
		warning("for efficiency of the bootstrap procedure, 'glm' should be called with 'x = TRUE'")
		break
	  }
	}
	
	rank <- match.arg(rank)
	
	ic <- switch(rank, AICc = AICc, AIC = AIC, BIC = BIC)
	
	mseq <- seq.int(M)
	
	if(all(vapply(models, inherits, FALSE, "glm"))) { ## !force.update &&
	  	best <- integer(R)
		ics <- numeric(M)
		no <- nobs(models[[1L]]) # assuming nobs is the same across models
		r <- 1L
		while(r <= R) {
			g <- sample.int(no, replace = TRUE)

			for(j in mseq) {
				fit <- models[[j]]
				fit1 <- glm.fit(model.matrix(fit)[g, , drop = FALSE], fit$y[g], fit$prior.weights[g],
				  offset = fit$offset[g], family = fit$family)
				ics[j] <- ic(loglik_glm_fit(fit1))
			}
			best[r] <- which.min(ics)
			r <- r + 1L
		}
	} else stop("all model objects must be of class \"glm\"")
  
	wts <- tabulate(best, M)
	wts <- wts  / sum(wts)
	structure(wts, wt.type = "bootstrap", names = names(models),
			  class = c("model.weights", class(wts)))
}


