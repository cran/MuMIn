`par.avg` <-
function(x, se, npar, weight, alpha = 0.05) {

	weight[is.na(weight)] <- 0 # not really necessary

	wx <- weighted.mean(x, weight, na.rm = TRUE)
	x.sqdiff <- (x - wx)^2
	xvar <- se^2

	avar <- weighted.mean(xvar + x.sqdiff, weight, na.rm = TRUE)^2
	ase <- weighted.mean(sqrt(xvar + x.sqdiff), weight, na.rm = TRUE)

	z <- c((qt(1 - (alpha / 2), npar) / qnorm(1 - (alpha / 2)))^2)

	use <- weighted.mean(sqrt((xvar * z) + x.sqdiff), weight, na.rm = TRUE)

	ci <- qnorm(1 - (alpha / 2)) * use

	return(c(`Coefficient` = wx, `Variance` = avar,  `SE` = ase, `Unconditional SE` = use, `Lower CI` = wx - ci, `Upper CI` = wx + ci))
}
