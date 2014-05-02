`coeffs` <-
function (model) UseMethod("coeffs")

`coeffs.gls` <-
function (model) summary(model)$coefficients

`coeffs.lme` <-
function(model) model$coefficients$fixed

# `coeffs.glmer` <-
# `coeffs.lmer` <-
# function(model) {
	# ret <- model@fixef
	# names(ret) <- model@cnames$.fixed
	# return(ret)
# }

`coeffs.mer` <-
function(model) model@fixef

`coeffs.merMod` <-
function (model) fixef(model)

`coeffs.coxme` <-
`coeffs.lmekin` <-
function(model) {
	# for class coxme:
	ret <- model$coefficients
	# for class lmekin and older coxme
	if(is.list(ret) && !is.null(ret$fixed)) return(ret$fixed)
	ret
}

`coeffs.unmarkedFit` <- 
function(model) {
	ret <- lapply(model@estimates@estimates, coef, altNames = FALSE)
	pfx <- rep(vapply(model@estimates@estimates, slot, "", "short.name"),
		vapply(ret, length, 1L))
	ret <- unlist(unname(ret))
	Ints <- which(names(ret) == "Int")
	names(ret) <- paste(pfx, "(", names(ret), ")", sep = "")
	attr(ret, "Intercept") <- Ints
	ret
}

`coeffs.splm` <- 
function (model) {
	c(model$coefficients, model$arcoef,
	if(is.matrix(model$errcomp)) model$errcomp[, 1L] else model$errcomp)
}

`coeffs.MCMCglmm` <-
function (model)
#summary(model)$solutions[, 1L]
colMeans(model$Sol[, seq.int(model$Fixed$nfl), drop = FALSE])



`coeffs.gamm` <-
function (model) coef(model$gam)

`coeffs.mark` <- function(model) {
	cf <- model$results$beta[, 1L]
	names(cf) <- gsub("^([a-zA-Z]+):(.*)$", "\\1(\\2)",
		rownames(model$results$beta), perl = TRUE)
	cf
}

`coeffs.multinom` <- 
function (model) {
	coefs <- coef(model)
	if (is.vector(coefs)) {
      	coefs <- t(as.matrix(coefs))
    	}
    	coefdim <- dim(coefs)
	Names <- dimnames(coefs)
	if (is.null(Names[[1L]])) 
      	Names <- Names[[2L]]
    	else Names <- as.vector(outer(Names[[2L]], Names[[1L]], function(name2, 
      	name1) paste(name1, name2, sep = ":")))
	res <- as.vector(coefs)
	names(res) <- Names
	res
}


`coeffs.aodml` <-
function (model) {
	c(model$b, model$phi)
}

`coeffs.asreml` <- 
function (model) {
	coef(model)$fixed  ## should include also '$sparse' ?
}

`coeffs.cpglmm` <-
function (model) 
model@fixef

`coeffs.default` <-
function(model) coef(model)
