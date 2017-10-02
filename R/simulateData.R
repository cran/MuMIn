.simulateData <-
function(n = 10, k = 2, family = gaussian,
		 beta = runif(k),
		 intercept = TRUE,
		 gamma.shape = 1, gaussian.sd = 1,
		 binomial.n = 1,
		  ...) {
	
	if(missing(k) && !missing(beta)) k <- length(beta)
	if (is.character(family)) family <- get(family, mode = "function",
		envir = parent.frame())
    if (is.function(family))  family <- family()
	nv <- k - if(intercept) 1L else 0L
	tmp <- numeric(k)
	tmp[] <- beta
	beta <- tmp
	
	vs <- seq.int(nv)
	dat <- matrix(ncol = nv, nrow = n)
	
	if(k > 1) {
		colnames(dat) <- sprintf("X%d", vs)
		for(i in vs) dat[, i] <- runif(n, min = 0, max = 1)
	}
	eta <- ((if(intercept) cbind(1, dat, deparse.level = 0L) else dat) %*% beta)[, 1L]
	mu <- family$linkinv(eta) #
	
	y <- switch(family$family,
		   gaussian = rnorm(n, mean = mu, sd = gaussian.sd),
		   binomial = {
			if(is.function(binomial.n)) binomial.n <- binomial.n(length(mu))
			binomial.n <- rep(binomial.n, length.out = length(mu))
			y <- rbinom(n, size = binomial.n, mu)
			if(all(binomial.n == 1)) y else
				cbind(y, binomial.n - y, deparse.level = 0)
			},
		   poisson = rpois(n, lambda = mu),
		   Gamma = rgamma(n, rate = gamma.shape / mu, shape = gamma.shape),
		   stop("'family' not recognized")
	)
	dat <- cbind(y = NA, as.data.frame(dat))
	dat$y <- y
	dat
}