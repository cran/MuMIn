# `asscount` stands, obviously, for "all-subsets count".

asscount <-
function(object) {
    .isssmat <- \(x) is.matrix(x) && {
            d <- dim(x)
            d[1L] == d[2L]
        } && is.logical(x) &&
        !any(x[lower.tri(x)])

    ssmat <- if(.isssmat(object)) object else
        if(is.character(object)) attr(object,"deps") else
            if(is.object(object)) attr(getAllTerms(object),"deps")

    if(!.isssmat(ssmat)) stop("cannot get subset matrix from 'object'")
    diag(ssmat) <- TRUE
    numcombin(ssmat)
}
 
numcombin <-
function(subsetmatrix) {
  
  storedresults <- new.env(parent = emptyenv())
  makekey <- function(idx) paste0(sort(idx), collapse = "_")
  
  # counter <- integer(2L)
  
  .countideals <-
  function(idx) {
    # counter[1L] <<- counter[1L] + 1L 

    if (length(idx) == 0L) return(1L)
    
    key <- makekey(idx)
    if (exists(key, storedresults, inherits = FALSE))
        return(storedresults[[key]])

    #counter[2L] <<- counter[2L] + 1L 
    
    submat <- subsetmatrix[idx, idx, drop = FALSE]
    
    # Find maximal elements correctly:
    # maximal = no STRICT superset in idx
    M.idx <- idx[match(TRUE, rowSums(submat & upper.tri(submat)) == 0L)]
    
    # Case 1: exclude M
    idx1 <- idx[idx != M.idx]
    # Case 2: include M -> remove all its subsets
    idx2 <- setdiff(idx, idx[subsetmatrix[idx, M.idx]])
    
    res <- Recall(idx1) + Recall(idx2)
    assign(key, res, envir = storedresults)
    return(res)
  }
  .countideals(seq_len(ncol(subsetmatrix)))
}
