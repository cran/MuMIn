
cos2Weights <- 
function(object, ..., data, eps = 1E-6, maxit = 100, predict.args = list()) {
	models <- getModelArgs()
	M <- length(models)
	if(M < 2) stop("need more than one model")
 
    if(!all(vapply(models, inherits, TRUE, "lm")))
	   stop("'models' must inherit from \"lm\" class")
  
  	#py <- sapply(models, predict, newdata = data, type = "response", ...)
	cl <- as.call(c(as.name("predict"), alist(models[[i]], newdata = data,
		type = "response"), predict.args))
	
	# TODO: glm.fit version
	
	i <- 1L
	py1 <- eval(cl)
	if(is.array(py1)) {
		stop(">1-dimensional predictions cannot be handled yet")
		# assuming prediction is a matrix: (add some checking for it)
		py <- array(dim = c(dim(py1), M))
		py[, , 1L] <- py1
		for(i in seq.int(2L, M)) py[, , i] <- eval(cl)
	} else {
		py <- array(dim = c(length(py1), M))
		py[, 1L] <- py1
		for(i in seq.int(2L, M)) py[, i] <- eval(cl)
	}
	
	# if one model is constant:
	if (any(g <- apply(py, 2L, "sd") == 0))
		py[, g] <- py[, g] + rnorm(NROW(py))
		
	sqrtm <- getFrom("expm", "sqrtm")
	
	R <- cor(py)
	nR <- NCOL(R)
	D1 <- diag(rep(2, nR))
	D2 <- diag(nR)
	counter <- 0L
	while (any(abs(diag(D1) - diag(D2)) > eps)) {
		ED <- eigen(D1 %*% R %*% D1)
		Q <- ED$vectors
		Lambda <- diag(ED$values)
		## test:
		#Q %*% Lambda %*% solve(Q) # fine
		Lambda12 <- sqrtm(Lambda)
		E <- solve(D1) %*% Q %*% Lambda12 %*% solve(Q)
		D2 <- D1
		D1 <- diag(diag(Re(E)))
		counter <- counter + 1L
		if (counter >= maxit) {
		  warning("maximum number of iterations reached without convergence")
		  break
		}
	}
	wts <- diag(D2)^2 / sum(diag(D2)^2)
	
	structure(wts, wt.type = "cos-squared", names = names(models),
			  class = c("model.weights", class(wts)))
}
