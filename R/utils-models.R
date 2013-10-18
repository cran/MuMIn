


# Hidden functions
# `.getLogLik` <- function() logLik

`.getLogLik` <- function()
	if(isGeneric("logLik")) .xget("stats4", "logLik") else
		.xget("stats", "logLik")
	#if(isGeneric("logLik")) stats4:::logLik else
	#	stats::logLik

		
`.getLik` <- function(x) {
    if(isGEE(x)) {
		logLik <- quasiLik
		lLName <- "qLik"
	} else {
	 	logLik <- .getLogLik()
		lLName <- "logLik"
    }
	list(logLik = logLik, name = lLName)
}

`.getRank` <- function(rank = NULL, rank.args = NULL, object = NULL, ...) {
	rank.args <- c(rank.args, list(...))

	if(is.null(rank)) {
		x <- NULL # just not to annoy R check
		IC <- as.function(c(alist(x =, do.call("AICc", list(x)))))
		attr(IC, "call") <- call("AICc", as.name("x"))
		class(IC) <- c("function", "ICWithCall")
		return(IC)
	} else if(inherits(rank, "ICWithCall") && length(rank.args) == 0L) {
		return(rank)
	}

	srank <- substitute(rank, parent.frame())
	if(srank == "rank") srank <- substitute(rank)

	rank <- match.fun(rank)
	ICName <- switch(mode(srank), call = as.name("IC"), character = as.name(srank), name=, srank)
	ICarg <- c(list(as.name("x")), rank.args)
	ICCall <- as.call(c(ICName, ICarg))
	IC <- as.function(c(alist(x =), list(substitute(do.call("rank", ICarg), 
		list(ICarg = ICarg)))))

	if(!is.null(object)) {
		test <- IC(object)
		if (!is.numeric(test) || length(test) != 1L)
			stop("'rank' should return numeric vector of length 1")
	}

	attr(IC, "call") <- ICCall
	class(IC) <- c("function", "ICWithCall")
	IC
}

`matchCoef` <- function(m1, m2,
	all.terms = getAllTerms(m2, intercept = TRUE),
	beta = FALSE,
	terms1 = getAllTerms(m1, intercept = TRUE),
	coef1 = if (beta) beta.weights(m1)[, 3L] else coeffs(m1),
	allCoef = FALSE,
	...
	) {
	if(any((terms1 %in% all.terms) == FALSE)) stop("'m1' is not nested within 'm2'")
	row <- structure(rep(NA_real_, length(all.terms)), names = all.terms)

	fxdCoefNames <- fixCoefNames(names(coef1))
	row[terms1] <- NaN
	pos <- match(terms1, fxdCoefNames, nomatch = 0L)
	row[fxdCoefNames[pos]] <- coef1[pos]
	if(allCoef) {
		ct <- coefTable(m1, ...)
		i <- match(names(coef1), rownames(ct))
		j <- !is.na(i)
		rownames(ct)[i[j]] <- fxdCoefNames[j]
		attr(row, "coefTable") <- ct
	}
	row
}


# sorts alphabetically interaction components in model term names
# if 'peel', tries to remove coefficients wrapped into function-like syntax
# (this is meant mainly for 'unmarkedFit' models with names such as "psi(a:b:c)")
# TODO: this function is ugly, do something with it.
`fixCoefNames` <- function(x, sort = FALSE, peel = TRUE) {
	if(!length(x)) return(x)
	ox <- x
	ia <- grep(":", x, fixed = TRUE)
	if(!length(ia)) return(structure(x, order = rep.int(1L, length(x))))
	x <- ret <- x[ia]
	if(peel) {
		# case of pscl::hurdle, cf are prefixed with count_|zero_
		if(all(substr(x, 1L, pos <- regexpr("_", x, fixed = TRUE)) %in%
			c("count_", "zero_"))) {
				ret <- substr(ret, pos + 1L, 256L)
				k <- TRUE
				suffix <- ""
		} else { # unmarkedFit with its phi(...), lambda(...) etc...
			k <- grepl("^\\w+\\(.+\\)$", x, perl = TRUE)
			fname <- substring(x[k], 1L, attr(regexpr("^\\w+(?=\\()", x[k],
				perl = TRUE),"match.length"))
			# exclude common transformations
			k[k] <- !(fname %in% c("log", "I", "exp", "s", "te"))
			if(any(k)) {
				pos <- vapply(x[k], function(z) {
					parens <- lapply(lapply(c("(", ")"),
						function(s) gregexpr(s, z, fixed = TRUE)[[1L]]),
							function(y) y[y > 0L])
					parseq <- unlist(parens, use.names = FALSE)
					p <- cumsum(rep(c(1L, -1L), sapply(parens, length))[order(parseq)])
					if(any(p[-length(p)] == 0L)) -1 else parseq[1L]
				}, numeric(1L), USE.NAMES = FALSE)
				k[k] <- pos != -1
				pos <- pos[pos != -1]
				if(any(k)) ret[k] <- substring(x[k], pos + 1L, nchar(x[k]) - 1L)
			}
			suffix <- ")"
		}
	} else	k <- FALSE
	ret <- vapply(lapply(spl <- strsplit(ret, ":"), base::sort), paste, "",
		collapse = ":")
	if(peel && any(k))
		ret[k] <- paste(substring(x[k], 1L, pos), ret[k], suffix, sep = "")
	ox[ia] <- ret
	ord <- rep.int(1, length(ox))
	ord[ia] <- sapply(spl, length)
	structure(ox, order = ord)
}

#Tries to find out whether the models are fitted to the same data
.checkModels <- function(models, error = TRUE) {
	#
	cl <- sys.call(sys.parent())
	err <-  if (error) 	function(x) stop(simpleError(x, cl))
		else function(x) warning(simpleWarning(x, cl))
	res <- TRUE

	responses <- lapply(models, function(x) {
	  f <- formula(x)
	  if((length(f) == 2L) || (is.call(f[[2L]]) && f[[2L]][[1L]] == "~")) 0 else f[[2L]]
	})

 	if(!all(vapply(responses[-1L], "==", logical(1), responses[[1L]]))) {
		err("response differs between models")
		res <- FALSE
	}

	datas <- lapply(models, function(x) .getCall(x)$data)
	# when using only 'nobs' - seems to be evaluated first outside of MuMIn namespace
	# which e.g. gives an error in glmmML - the glmmML::nobs method is faulty.
	nresid <- vapply(models, function(x) nobs(x), numeric(1L)) # , nall=TRUE

	if(!all(datas[-1L] == datas[[1L]]) || !all(nresid[-1L] == nresid[[1L]])) {
		err("models are not all fitted to the same data")
		res <- FALSE
	}
	invisible(res)
}

#system.time(for(i in 1:1000) abbreviateTerms(x))

`abbreviateTerms` <- function(x, minlength = 4, minwordlen = 1,
	capwords = FALSE, deflate = FALSE) {
	if(!length(x)) return(x)
	if(deflate) dx <-
		#gsub("([\\(,]) *\\w+ *= *(~ *(1 *[\\+\\|]?)?)? *", "\\1", x, perl = TRUE)
		gsub("([,\\(\\[]|^)( *~ *)(1 *([\\|\\+] *)?)?", "\\1",
			gsub("([\\(,]) *\\w+ *= *", "\\1", x, perl = TRUE), perl = TRUE)
		else dx <- x

	#DebugPrint(x)
	s <- strsplit(dx, "(?=[\\W_])", perl = TRUE)
	# remove I(...):
	s <- lapply(s, function(z) {
		z <- if((n <- length(z)) > 3L && all(z[c(1L, 2L, n)] == c("I", "(", ")")))
			z[3L:(n - 1L)] else z
		z[z != " "]
	})
	v <- unique(unlist(s, use.names = FALSE))
	i <- grep("[[:alpha:]]", v, perl = FALSE)
	av <- v
	if(length(i)) {
		tb <- rbindDataFrameList(lapply(s, function(x)
			as.data.frame(rbind(c(table(x))))))
		tb[is.na(tb)] <- 0L

		if(length(v) > length(i)) minlength <-
			minlength - max(c(0L, apply(tb[, v[-i], drop = FALSE], 1L,
			"*", nchar(colnames(tb[, v[-i], drop = FALSE])))))
		n <- min(minlength / rowSums(tb[, v[i], drop = FALSE]))
		if(deflate) {
			repl1 <- c("TRUE" = "T", "FALSE" = "F", "NULL" = "")
			for(j in seq_along(repl1)) av[av == names(repl1)[j]] <- repl1[j]
		}
		av[i] <- abbreviate(av[i], max(n, minwordlen))
		if(capwords) av[i] <- paste(toupper(substring(av[i], 1L, 1L)),
				tolower(substring(av[i], 2L)), sep = "")
	}
	for(j in seq_along(s)) s[[j]] <- paste(av[match(s[[j]], v)], collapse = "")
	names(av) <- v
	structure(unlist(s), names = x, variables = av[i])
}

`model.names` <- function(object, ..., labels = NULL, use.letters = FALSE) {
	if (missing(object) && length(models <- list(...)) > 0L) {
		object <- models[[1L]]
	} else if (inherits(object, "list")) {
		if(length(object) ==  0L) stop("at least one model must be given")
		models <- object
		object <- models[[1L]]
	} else models <- list(object, ...)
	if(length(models) == 0L) stop("at least one model must be given")
	.modelNames(models = models, uqTerms = labels, use.letters = use.letters)
}

`.modelNames` <- function(models = NULL, allTerms, uqTerms, use.letters = FALSE, ...) {
	if(missing(allTerms)) allTerms <- lapply(models, getAllTerms)
	if(missing(uqTerms) || is.null(uqTerms))
		uqTerms <- unique(unlist(allTerms, use.names = FALSE))
	
	n <- length(uqTerms)
	
	if(use.letters && n > length(LETTERS)) stop("more terms than there are letters")
	sep <- if(!use.letters && n > 9L) "/" else ""
	
	labels <- if (use.letters) LETTERS[seq_len(n)] else as.character(seq_len(n))
	ret <- sapply(allTerms, function(x) paste(labels[sort(match(x, uqTerms))],
		collapse = sep))

	dup <- table(ret)
	dup <- dup[dup > 1L]

	if(length(dup) > 0L) {
		idup <- which(ret %in% names(dup))
		ret[idup] <- sapply(idup, function(i) paste(ret[i],
			letters[sum(ret[seq.int(i)] == ret[i])], sep = ""))
	}
	ret[ret == ""] <- "(Null)"
	attr(ret, "variables") <- structure(seq_along(uqTerms), names = uqTerms)
	ret
}

`modelDescr` <- function(models, withModel = FALSE, withFamily = TRUE,
	withArguments = TRUE, remove.cols = c("formula", "random", "fixed", "model",
	"data", "family", "cluster", "model.parameters")) {

	if(withModel) {
		allTermsList <- lapply(models, function(x) {
			tt <- getAllTerms(x)
			rtt <- attr(tt, "random.terms")
			c(tt, if(!is.null(rtt)) paste("(", rtt, ")", sep = "") else NULL)
		})
		allTerms <- unique(unlist(allTermsList))
		abvtt <- abbreviateTerms(allTerms)
		variables <- attr(abvtt, "variables")
		abvtt <- gsub("\\(1 \\| (\\S+)(?: %in%.*)?\\)", "(\\1)", abvtt, perl = TRUE)
		abvtt <- sapply(allTermsList, function(x) paste(abvtt[match(x, allTerms)],
			collapse = "+"))
	} else abvtt <- variables <- NULL


	if(withFamily) {
		fam <- sapply(models, function(x) tryCatch(unlist(family(x)[c("family",
			"link")]), error = function(e) character(2L)) )

		f <- fam[1L, ]
		f[is.na(f)] <- ""
		f <- vapply(strsplit(f, "(", fixed = TRUE), "[", "", 1L)
		f[f == "Negative Binomial"] <- "negative.binomial"

		fam[2L, fam[2L, ] == vapply(unique(f), function(x) if(is.na(x))
									NA_character_ else formals(get(x))$link,
			FUN.VALUE = "")[f]] <- NA_character_

		j <- !is.na(fam[2L,])
		fnm <- fam[1L, j]
		fnm <- ifelse(substring(fnm, nchar(fnm)) != ")",
			paste(fnm, "(", sep = ""), paste(substring(fnm, 1, nchar(fnm) - 1),
				", ", sep = ""))
		fam[1L, j] <- paste(fnm, fam[2L, j], ")", sep = "")
	}

	if(withArguments) {
		cl <- lapply(models, .getCall)
		haveNoCall <-  vapply(cl, is.null, FALSE)
		cl[haveNoCall] <- lapply(cl[haveNoCall], function(x) call("none", formula = NA))
 		arg <- lapply(cl, function(x) sapply(x[-1L], function(argval)
			switch(mode(argval), character = , logical = argval,
			numeric = signif(argval, 3L), deparse(argval, nlines = 1L))))
		arg <- rbindDataFrameList(lapply(lapply(arg, t), as.data.frame))
		arg <- cbind(class = as.factor(sapply(lapply(models, class), "[", 1L)),
			arg[, !(colnames(arg) %in% remove.cols), drop = FALSE])
		reml <-	rep(NA, length(models))
		if(!is.null(arg$method)) {
			reml <- ((arg$class == "lme" &
				is.na(arg$method)) | arg$method == "REML")
			arg$method  <- NULL
		}
		if(!is.null(arg$REML)) reml <- ifelse(is.na(arg$REML), reml, arg$REML == "TRUE")
		arg$REML <- as.factor(reml)

		arg <- as.matrix(arg)
		arg[is.na(arg) | arg == "NULL"] <- ""
		arg <- arg[, apply(arg, 2L, function(x) length(unique(x))) != 1L, drop = FALSE]
		if(ncol(arg)) arg <- gsub("([\"'\\s]+|\\w+ *=)","", arg, perl = TRUE)
	}
	ret <- as.data.frame(cbind(model = abvtt, family = fam[1L, ], arg, deparse.level = 0L))
	attr(ret, "variables") <- variables
	ret
}
