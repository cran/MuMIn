`par.avg` <-
function(x, se, weight, df = NULL, alpha = 0.05, revised.var = TRUE) {

	#print(match.call())

	if (!(is.numeric(x) && is.numeric(se) && is.numeric(weight)))
		stop("'x', 'se' and 'weight' must be numeric vectors")
	n <- length(x)
	if(length(weight) != n || length(se) != n)
		stop("'x', 'se' and 'weight' are not of the same length: ",
			 sprintf("x: %d, weight: %d, se: %d", n, length(weight),
			 length(se)))

	weight[is.na(weight)] <- 0 # not really necessary

	wx <- weighted.mean(x, weight, na.rm = TRUE)
	x.sqdiff <- (x - wx)^2
	xvar <- se^2

	do.ase <- !(missing(df) || is.null(df) || all(is.na(df)))

	if(do.ase) {
		z <- c((qt(1 - (alpha / 2), df) / qnorm(1 - (alpha / 2)))^2)
		#print(cbind(x, z, weight))

		i <- is.na(df) & !is.na(x)
		if (length(i) > 0) z[i] <- 0

		#print(cbind(x, z, weight))
	}
	#print(rbind(xvar, weight))

	if(revised.var) {
		# Unconditional sqrt-Variance, revised in B&A2004
		use <- sqrt(weighted.mean(xvar + x.sqdiff, weight, na.rm = TRUE))
		if (do.ase)
			# Adjusted std. error - formula modified analogically to the previous
			ase <- sqrt(weighted.mean((xvar * z) + x.sqdiff, weight, na.rm = TRUE))
	} else {
		# Unconditional sqrt-Variance, original formula (B&A2002, eqn 4.7)
		use <- weighted.mean(sqrt(xvar + x.sqdiff), weight, na.rm = TRUE)
		if (do.ase)
			# Adjusted std. error (B&A2002, ch.4.3.3/p164):
			ase <- weighted.mean(sqrt((xvar * z) + x.sqdiff), weight, na.rm = TRUE)
	}

	ci <- qnorm(1 - (alpha / 2)) * (if (do.ase) ase else use)

	return(c(`Coefficient` = wx, `SE` = use,
			 `Adjusted SE` = if(do.ase) ase else NA,
			 `Lower CI` = wx - ci, `Upper CI` = wx + ci))
}
