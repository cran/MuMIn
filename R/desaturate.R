
desaturate <-
function(x, f = 1) {
    r <- x[, 1L]
    g <- x[, 2L]
    b <- x[, 3L]
    L <- 0.3 * r + 0.6 * g + 0.1 * b
    #L <- 0.6 * r + 0.1 * g + 0.3 * b
    cbind(r + f * (L - r),
          g + f * (L - g),
          b + f * (L - b))
}
