`formulaAllowed` <-
function(frm) {
	factors <- attr(terms(frm), "factors")
	return(all(factors < 2))
}

