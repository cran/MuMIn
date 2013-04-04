`plot.model.selection` <-
function(x, col = c("SlateGray", "SlateGray3"),
	ylab=expression("Cumulative Akaike weight" ~~ (omega)),
	xlab=NA, col.bias=3, col2="white", ...) {
	cumw <- cumsum(x$weight)
	n <- nrow(x)
	v <- attr(x, "terms")
	m <- length(v)
	plot(c(0, m), c(0, 1), type="n", axes=F, ann=F,
		ylim = c(1, 0),...)
	#col <- rep(col, length.out=m)
	#colp <- colorRamp(c("SlateGrey","black"))

	pal <- lapply(col, function(x) colorRamp(c(col2, x, x), bias=col.bias))
	npal <- length(col)

	for(i in 1:m) {
		rect(i - 1, c(0, cumw), i, c(cumw, 1),
			#col=ifelse(is.na(x[,i]), "white", col[i])
			col=ifelse(is.na(x[,i]), "white",
			rgb(pal[[1L + ((i - 1L) %% npal)]](x$weight)/255))
		)
	}
	#mtext(side=3, at=1:m - 0.5, v)
	lapply(1:m, function(i)	mtext(parse(text=v[[i]]), side=3, at=i- 0.5,
		padj=0.5, line=1))

	axis(2)
	ss <- x$weight > -strheight("I")
	mtext(side=4, at=(c(0, cumw[-n]) + cumw)[ss]/2, rownames(x)[ss],
		las=1, line=1.25, adj=1)
	title(ylab=ylab, xlab=xlab)
}

#plot(dd12, col = c("SlateGray", "SlateGray3"))
