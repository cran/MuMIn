umf_specs <-
structure(list(unmarkedFitColExt = list(structure(c("psiformula", 
"psi", "psi", "gammaformula", "col", "col", "epsilonformula", 
"ext", "ext", "pformula", "p", "det"), .Dim = 3:4, .Dimnames = list(
    c("argument", "est:short.name", "est:label"), c("psi", "gam", 
    "eps", "p")), umf_class = "unmarkedFitColExt", n_formulas = 4)), 
    unmarkedFitDS = list(structure(c("formula", "p", "det", "formula", 
    "lam", "state"), .Dim = c(3, 2), .Dimnames = list(c("argument", 
    "est:short.name", "est:label"), c("p", "lam")), umf_class = "unmarkedFitDS", n_formulas = 2, reversed = TRUE), 
        structure(c("formula", "", "", "formula", "lam", "state"
        ), .Dim = c(3, 2), .Dimnames = list(c("argument", "est:short.name", 
        "est:label"), c("dummy1", "lam")), umf_class = "unmarkedFitDS", n_formulas = 2, reversed = TRUE), 
        structure(c("formula", "lam", "state", "formula", "p", 
        "scale"), .Dim = c(3, 2), .Dimnames = list(c("argument", 
        "est:short.name", "est:label"), c("lam", "p")), umf_class = "unmarkedFitDS", n_formulas = 2)), 
    unmarkedFitGDS = list(structure(c("lambdaformula", "lambda", 
    "lambda", "phiformula", "phi", "phi", "pformula", "p", "det"
    ), .Dim = c(3, 3), .Dimnames = list(c("argument", "est:short.name", 
    "est:label"), c("lam", "phi", "p")), umf_class = "unmarkedFitGDS", n_formulas = 3), 
        structure(c("lambdaformula", "lambda", "lambda", "phiformula", 
        "phi", "phi", "pformula", "p", "det", "", "scale", "scale"
        ), .Dim = 3:4, .Dimnames = list(c("argument", "est:short.name", 
        "est:label"), c("lam", "phi", "p", "scale")), umf_class = "unmarkedFitGDS", n_formulas = 3), 
        structure(c("lambdaformula", "lambda", "lambda", "phiformula", 
        "phi", "phi", "pformula", "", ""), .Dim = c(3, 3), .Dimnames = list(
            c("argument", "est:short.name", "est:label"), c("lam", 
            "phi", "p")), umf_class = "unmarkedFitGDS", n_formulas = 3)), 
    unmarkedFitGMM = list(structure(c("lambdaformula", "lambda", 
    "lambda", "phiformula", "phi", "phi", "pformula", "p", "det"
    ), .Dim = c(3, 3), .Dimnames = list(c("argument", "est:short.name", 
    "est:label"), c("lam", "phi", "p")), umf_class = "unmarkedFitGMM", n_formulas = 3)), 
    unmarkedFitGPC = list(structure(c("lambdaformula", "lambda", 
    "lambda", "phiformula", "phi", "phi", "pformula", "p", "det"
    ), .Dim = c(3, 3), .Dimnames = list(c("argument", "est:short.name", 
    "est:label"), c("lam", "phi", "p")), umf_class = "unmarkedFitGPC", n_formulas = 3)), 
    unmarkedFitMPois = list(structure(c("formula", "p", "det", 
    "formula", "lambda", "state"), .Dim = c(3, 2), .Dimnames = list(
        c("argument", "est:short.name", "est:label"), c("p", 
        "lam")), umf_class = "unmarkedFitMPois", n_formulas = 2, reversed = TRUE)), 
    unmarkedFitOccu = list(structure(c("formula", "p", "det", 
    "formula", "psi", "state"), .Dim = c(3, 2), .Dimnames = list(
        c("argument", "est:short.name", "est:label"), c("p", 
        "psi")), umf_class = "unmarkedFitOccu", n_formulas = 2, reversed = TRUE)), 
    unmarkedFitOccuPEN = list(structure(c("formula", "p", "det", 
    "formula", "psi", "state"), .Dim = c(3, 2), .Dimnames = list(
        c("argument", "est:short.name", "est:label"), c("p", 
        "psi")), umf_class = "unmarkedFitOccuPEN", n_formulas = 2, reversed = TRUE)), 
    unmarkedFitOccuPEN_CV = list(structure(c("formula", "p", 
    "det", "formula", "psi", "state"), .Dim = c(3, 2), .Dimnames = list(
        c("argument", "est:short.name", "est:label"), c("p", 
        "psi")), umf_class = "unmarkedFitOccuPEN_CV", n_formulas = 2, reversed = TRUE)), 
    unmarkedFitOccuRN = list(structure(c("formula", "p", "det", 
    "formula", "lam", "state"), .Dim = c(3, 2), .Dimnames = list(
        c("argument", "est:short.name", "est:label"), c("p", 
        "lam")), umf_class = "unmarkedFitOccuRN", n_formulas = 2, reversed = TRUE)), 
    unmarkedFitPCO = list(structure(c("lambdaformula", "lam", 
    "lambda", "gammaformula", "gam", "gamma", "omegaformula", 
    "omega", "omega", "pformula", "p", "det", "iotaformula", 
    "", ""), .Dim = c(3, 5), .Dimnames = list(c("argument", "est:short.name", 
    "est:label"), c("lam", "gam", "om", "p", "iota")), umf_class = "unmarkedFitPCO", n_formulas = 5), 
        structure(c("lambdaformula", "lam", "lambda", "gammaformula", 
        "", "", "omegaformula", "omega", "omega", "pformula", 
        "p", "det", "iotaformula", "", ""), .Dim = c(3, 5), .Dimnames = list(
            c("argument", "est:short.name", "est:label"), c("lam", 
            "gam", "om", "p", "iota")), umf_class = "unmarkedFitPCO", n_formulas = 5), 
        structure(c("lambdaformula", "lam", "lambda", "gammaformula", 
        "gam", "gamma", "omegaformula", "", "", "pformula", "p", 
        "det", "iotaformula", "", ""), .Dim = c(3, 5), .Dimnames = list(
            c("argument", "est:short.name", "est:label"), c("lam", 
            "gam", "om", "p", "iota")), umf_class = "unmarkedFitPCO", n_formulas = 5), 
        structure(c("lambdaformula", "lam", "lambda", "gammaformula", 
        "gam", "gamma", "omegaformula", "om", "omega", "pformula", 
        "p", "det", "iotaformula", "", ""), .Dim = c(3, 5), .Dimnames = list(
            c("argument", "est:short.name", "est:label"), c("lam", 
            "gam", "om", "p", "iota")), umf_class = "unmarkedFitPCO", n_formulas = 5)), 
    unmarkedFitPCount = list(structure(c("formula", "p", "det", 
    "formula", "lam", "state", "", "alpha", "alpha"), .Dim = c(3, 
    3), .Dimnames = list(c("argument", "est:short.name", "est:label"
    ), c("p", "lam", "alpha")), umf_class = "unmarkedFitPCount", n_formulas = 2, reversed = TRUE), 
        structure(c("formula", "p", "det", "formula", "lam", 
        "state"), .Dim = c(3, 2), .Dimnames = list(c("argument", 
        "est:short.name", "est:label"), c("p", "lam")), umf_class = "unmarkedFitPCount", n_formulas = 2, reversed = TRUE))), .Names = c("unmarkedFitColExt", 
"unmarkedFitDS", "unmarkedFitGDS", "unmarkedFitGMM", "unmarkedFitGPC", 
"unmarkedFitMPois", "unmarkedFitOccu", "unmarkedFitOccuPEN", 
"unmarkedFitOccuPEN_CV", "unmarkedFitOccuRN", "unmarkedFitPCO", 
"unmarkedFitPCount"))
