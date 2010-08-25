`dredge` <-
function(global.model, beta = FALSE, eval = TRUE, rank = "AICc",
		 fixed = NULL, m.max = NA, subset, ...) {

	rankFn <- match.fun(rank)
	if (is.function(rank))
  		rank <- deparse(substitute(rank))

	rank.custom <- !missing(rank)
	if (rank.custom) {
		arg <- list(...)
		rankFnCall <- as.call(c(as.name("rankFn"), as.symbol("x"), arg))
		IC <- function(x) eval(rankFnCall)
		res <- IC(global.model)
  		if (!is.numeric(res) || length(res) != 1) {
			stop("'rank' should return numeric vector of length 1")
		}
	} else {
		rankFnCall <- call("AICc", as.symbol("x"))
	}
	#print(rankFnCall)

	intercept <- "(Intercept)"

	all.terms <- getAllTerms(global.model)

	has.int <- attr(all.terms, "intercept")
	tmp <- attributes(all.terms)
	all.terms <- fixCoefNames(all.terms)
	attributes(all.terms) <- tmp

	n.vars <- length(all.terms)
	ms.tbl <- numeric(0)
	formulas <- character(0)

	is.glm <- inherits(global.model, "glm")
	is.lm <- !is.glm & inherits(global.model, "lm")

	if (
			(inherits(global.model, c("mer")) && ("REML" %in%
				names(deviance(global.model))))
		|| 	(inherits(global.model, c("lme")) && global.model$method == "REML")
		||  (any(inherits(global.model, c("lmer", "glmer")))
			  && global.model@status["REML"] != 0)
	) {
		warning("Comparing models with different fixed effects fitted by REML")
	}

	if (beta && is.null(tryCatch(beta.weights(global.model), error=function(e) NULL,
		warning=function(e) NULL))) {
		warning("Do not know how to calculate B-weigths for ",
				class(global.model)[1], ", argument 'beta' ignored")
         beta <- FALSE
	}

	summary.globmod <- summary(global.model)

	has.rsq <- is.list(summary.globmod) && is.numeric(summary.globmod$r.squared)
	has.dev <- !is.null(deviance(global.model))

	m.max <- if (missing(m.max)) n.vars else min(n.vars, m.max)

	# fixed variables:
	if (!is.null(fixed)) {
		if (inherits(fixed, "formula")) {
			if (fixed[[1]] != "~" || length(fixed) != 2)
				warning("'fixed' attribute should be a formula of form: ",
						"~ a + b + c")
			fixed <- c(getAllTerms(fixed))
		} else if (!is.character(fixed)) {
			stop ("'fixed' should be either a character vector with"
				  + " names of variables or a one-sided formula")
		}
		if (!all(fixed %in% all.terms)) {
			warning("Not all terms in 'fixed' exist in 'global.model'")
			fixed <- fixed[fixed %in% all.terms]
		}
		m.max <- m.max - length(fixed)
	}

	if (m.max > 0) {
		num.opt.vars <- 1:n.vars
		if (!is.null(fixed))
			num.opt.vars <- num.opt.vars[!(all.terms %in% fixed)]

		all.comb <- lapply(seq(m.max), combn, x = num.opt.vars)
		all.comb <- unlist(lapply(all.comb, function(.x) split(.x, col(.x))),
						   recursive = FALSE)
		all.comb <- c(`0` = list(0), all.comb)

	} else {
		all.comb <- list(0)
	}

	if (!is.null(fixed))
		all.comb <- lapply(all.comb, append, (1:n.vars)[all.terms %in% fixed])

	formulas <- lapply(all.comb, function(.x) reformulate(c("1", all.terms[.x]),
														  response = "." ))

	formulas <- lapply(formulas, `attr<-`, ".Environment",
				 attr(formula(global.model), ".Environment"))

	ss <- sapply(formulas, formulaAllowed)

	all.comb <- all.comb[ss]
	formulas <- formulas[ss]

	## Apply subset:
	if (!missing(subset)) {
		xtable <- as.data.frame(!is.na(do.call("rbind", lapply(all.comb, match,
			x = seq.int(n.vars)))))
		colnames(xtable) <- all.terms
		rownames(xtable) <- NULL

		ss <- eval(substitute(subset), envir = xtable)
		formulas <- formulas[ss]
		all.comb <- all.comb[ss]
		# 10 time slower! #xtable <- !is.na(t(sapply(all.comb, match, x=seq.int(n.vars))))
	}

	names(formulas) <- seq(formulas)
	if (any(inherits(global.model, c("mer", "lmer", "glmer")))) {
          formulas <- lapply(formulas, update, attr(all.terms, "random"))
	}

	if (!eval) return(formulas)

	getK <- function(x) as.vector(attr(logLik(x), "df"))

	### BEGIN:
	for(b in seq(length(all.comb))) {
		# print(all.comb[[b]])
        terms1 <- all.terms[all.comb[[b]]]
		frm <- formulas[[b]]

		row1 <- rep(NA, n.vars)
		row1[match(terms1, all.terms)] <- rep(1, length(terms1))

		cl <- call("update", substitute(global.model), frm)

		#TODO: optional error printing.
		fit1 <- tryCatch(eval(cl, parent.frame()), error=function(err) {
			err$message <- paste(conditionMessage(err), "(model skipped by dredge)")
			warning(err)
			return(NULL)
		})

		if (is.null(fit1)) {
			formulas[[as.character(b)]] <- NA
			next;
		}

		coef1 <- if (beta) beta.weights(fit1)[,3] else coeffs(fit1)
		names(coef1) <- fixCoefNames(names(coef1))

		icept <- if (has.int) coef1[intercept] else NA

		coef1 <- coef1[c(na.omit(match(all.terms, names(coef1))))]
		row1[match(names(coef1), all.terms)] <- coef1

		row1 <- c(icept, row1, k=getK(fit1))
		if (has.rsq) {
			fit1.summary <- summary(fit1)
			row1 <- c(row1, r.squared=fit1.summary$r.squared,
				adj.r.squared=fit1.summary$adj.r.squared)
		}
		if (has.dev)
			row1 <- c(row1, deviance(fit1))

		if (rank.custom) {
			ic <- IC(fit1)
			row1 <- c(row1, IC=ic)
		} else {
			aicc <- AICc(fit1)
		    row1 <- c(row1, AIC=attr(aicc, "AIC"), AICc=aicc)
		}
		ms.tbl <- rbind(ms.tbl, row1)
	} ### END

	formulas[is.na(formulas)] <- NULL
	ms.tbl <- data.frame(ms.tbl, row.names=seq(NROW(ms.tbl)))

	# Convert columns with presence/absence of terms to factors
	tfac <- which(c(FALSE, !(all.terms %in% fixCoefNames(names(coeffs(global.model))))))
	ms.tbl[tfac] <- lapply(ms.tbl[tfac], factor, levels=1, labels="+")


	colnames(ms.tbl) <- c("(int.)", all.terms, "k",
		if (has.rsq) c("R.sq", "Adj.R.sq"),
		if (has.dev) if (is.lm) "RSS" else "Dev.",
		if (rank.custom) rank else c("AIC", "AICc")
	)

	o <- order(ms.tbl[, rank], decreasing = FALSE)
	ms.tbl <- ms.tbl[o,]
	ms.tbl$delta <- ms.tbl[, rank] - min(ms.tbl[, rank])
	ms.tbl$weight <- exp(-ms.tbl$delta / 2) / sum(exp(-ms.tbl$delta / 2))

	class(ms.tbl) = c("model.selection", "data.frame")

	attr(ms.tbl, "formulas") <- formulas[o]
	attr(ms.tbl, "global") <- global.model
	attr(ms.tbl, "terms") <- c(intercept, all.terms)

	#print(rankFnCall)

	if (rank.custom) {
		rankFnCall[[1]] <- as.name(rank)
	}
	#rankFnCall[[2]] <- substitute(global.model)
	attr(ms.tbl, "rank.call") <- rankFnCall


	if (!is.null(attr(all.terms, "random.terms")))
		attr(ms.tbl, "random.terms") <- attr(all.terms, "random.terms")

	return(ms.tbl)
}

`subset.model.selection` <-
function(x, subset, select, recalc.weights = TRUE, ...) {
	if (missing(select)) {
		e <- substitute(subset)
		i <- eval(e, x, parent.frame(2))
		return(`[.model.selection`(x, i, recalc.weights = recalc.weights, ...))
	} else {
		cl <- match.call(expand.dots = FALSE)
	    cl <- cl[c(1L, match(names(formals("subset.data.frame")), names(cl), 0L))]
	    cl[[1L]] <- as.name("subset.data.frame")
	    return(eval(cl, parent.frame(2)))
	}
}

`[.model.selection` <-
function (x, i, j, recalc.weights = TRUE, ...) {
	res <- `[.data.frame`(x, i, j, ...)
	if (missing(j)) {
		attr(res, "global") <- attr(x, "global")
		attr(res, "terms") <- attr(x, "terms")
		attr(res, "formulas") <- attr(x, "formulas")[i]
		if(recalc.weights)
			res$weight <- res$weight / sum(res$weight)
		class(res) <- class(x)
	} else {
		class(res) <- "data.frame"
	}
	return(res)
}

`print.model.selection` <-
function(x, abbrev.names = TRUE, ...) {
	if(!is.null(x$weight))
		x$weight <- round(x$weight, 3)

	xterms <- attr(x, "terms")

	if(is.null(xterms)) {
		print.data.frame(x, ...)
	} else {

		xterms <- gsub(" ", "", xterms)

		if(abbrev.names) {
			xterms <- sapply(xterms, function(z) {
				spl <- strsplit(rep(z, 2), c("[^\\w\\.]+",  "[\\w\\.$]+"),
					perl=TRUE)
				spl[[1]] <- abbreviate(spl[[1]], minlength=3)
				o <- order(sapply(lapply(spl, `==`, ""), `[[`, 1), decreasing=TRUE)
				x2 <- unsplit(spl, rep(order(sapply(lapply(spl, `==`, ""),
					`[[`, 1), decreasing=T), length.out=length(unlist(spl))))
				return(paste(x2, collapse=""))
			})
		}

		colnames(x)[seq_along(xterms)] <-  xterms


		gm <- attr(x, "global")
		gm.call <- (if(mode(gm) == "S4") `@` else `$`)(gm, "call")
		if(!is.null(call)) {
			cat("Global model: ")
			print(gm.call)
		}

		cat ("---\nModel selection table \n")
		i <- sapply(x, is.numeric)
		x[,i] <- signif(x[,i], 4)
		print.default(as.matrix(x), na.print="", quote=FALSE)
		if (!is.null(attr(x, "random.terms"))) {
			cat("Random terms:", paste(attr(x, "random.terms"), collapse=", "),
				"\n")
		}
	}
}


#sorts alphabetically interaction components in model term names
`fixCoefNames` <-
function(x) {
	if(!is.character(x)) return(x)
	return(sapply(lapply(strsplit(x, ":"), sort), paste, collapse=":"))
}
