.modelspecs <-
list(`<default>` = structure(list(model = "mean", item.name = NA_character_, 
    formula.arg = "formula"), row.names = c(NA, -1L), class = c("model.specs", 
"data.frame"), object.info = c(package = "<default>", func = "<default>", 
className = "<default>")), betareg = structure(list(model = c("mean", 
"dispersion"), item.name = c("mean", "precision"), formula.arg = c("formula1", 
"formula2")), row.names = c(NA, -2L), class = c("model.specs", 
"data.frame"), object.info = c(package = "betareg", func = "betareg", 
className = "betareg")), gamlss = structure(list(model = c("mean", 
"dispersion", "power", "tau"), item.name = c("mu", "sigma", "nu", 
"tau"), formula.arg = c("formula", "sigma.formula", "nu.formula", 
"tau.formula")), row.names = c(NA, -4L), class = c("model.specs", 
"data.frame"), object.info = c(package = "gamlss", func = "gamlss", 
className = "gamlss")), glmmTMB = structure(list(model = c("mean", 
"zeroinfl", "dispersion"), item.name = c("cond", "zi", "disp"
), formula.arg = c("formula", "ziformula", "dispformula")), row.names = c(NA, 
-3L), class = c("model.specs", "data.frame"), object.info = c(package = "glmmTMB", 
func = "glmmTMB", className = "glmmTMB")), glmerMod = structure(list(
    model = "mean", item.name = NA_character_, formula.arg = "formula"), row.names = c(NA, 
-1L), class = c("model.specs", "data.frame"), object.info = c(package = "lme4", 
func = "glmer", className = "glmerMod")), lmerMod = structure(list(
    model = "mean", item.name = NA_character_, formula.arg = "formula"), row.names = c(NA, 
-1L), class = c("model.specs", "data.frame"), object.info = c(package = "lme4", 
func = "lmer", className = "lmerMod")), gls = structure(list(
    model = "mean", item.name = NA_character_, formula.arg = "model"), row.names = c(NA, 
-1L), class = c("model.specs", "data.frame"), object.info = c(package = "nlme", 
func = "gls", className = "gls")), lme = structure(list(model = "mean", 
    item.name = NA_character_, formula.arg = "fixed,random"), row.names = c(NA, 
-1L), class = c("model.specs", "data.frame"), object.info = c(package = "nlme", 
func = "lme", className = "lme")), nlme = structure(list(model = "mean", 
    item.name = NA_character_, formula.arg = "model,fixed,random"), row.names = c(NA, 
-1L), class = c("model.specs", "data.frame"), object.info = c(package = "nlme", 
func = "nlme", className = "nlme")), hurdle = structure(list(
    model = c("mean", "zeroinfl"), item.name = c("count", "zero"
    ), formula.arg = c("formula1", "formula2")), row.names = c(NA, 
-2L), class = c("model.specs", "data.frame"), object.info = c(package = "pscl", 
func = "hurdle", className = "hurdle")), zeroinfl = structure(list(
    model = c("mean", "zeroinfl"), item.name = c("count", "zero"
    ), formula.arg = c("formula1", "formula2")), row.names = c(NA, 
-2L), class = c("model.specs", "data.frame"), object.info = c(package = "pscl", 
func = "zeroinfl", className = "zeroinfl")), lm = structure(list(
    model = "mean", item.name = NA_character_, formula.arg = "formula"), row.names = c(NA, 
-1L), class = c("model.specs", "data.frame"), object.info = c(package = "stats", 
func = "lm", className = "lm")), unmarked = list(unmarkedFitColExt = list(
    colext = structure(list(model = c("initial occupancy", "colonization", 
    "extinction", "detection"), item.name = c("psi", "col", "ext", 
    "det"), formula.arg = c("psiformula", "gammaformula", "epsilonformula", 
    "pformula"), name = c("Initial", "Colonization", "Extinction", 
    "Detection"), short.name = c("psi", "col", "ext", "p"), formula.slot = c("psiformula", 
    "gamformula", "epsformula", "detformula"), formlist.name = c(NA_character_, 
    NA_character_, NA_character_, NA_character_), formula.slot2 = c("formula1", 
    "formula2", "formula3", "formula4"), data = structure(c(1L, 
    2L, 2L, 3L), levels = c("site", "site+year", "site+year+obs"
    ), class = "factor")), row.names = c(NA, -4L), class = c("model.specs", 
    "data.frame"), object.info = c(package = "unmarked", func = "colext", 
    className = "unmarkedFitColExt", fitType = "colext"))), unmarkedFitDS = list(
    distsamp = structure(list(model = c("abundance", "detection", 
    "scale"), item.name = c("state", "det", "scale"), formula.arg = c("formula2", 
    "formula1", NA), name = c("Abundance|Density", "Detection", 
    "Hazard-rate(scale)"), short.name = c("lam", "p", "p"), formula.slot = c("formula2", 
    "formula1", NA), formlist.name = c(NA_character_, NA_character_, 
    NA_character_), formula.slot2 = c(NA_character_, NA_character_, 
    NA_character_), data = structure(c(1L, 2L, NA), levels = c("site", 
    "site+obs"), class = "factor")), row.names = c(NA, -3L), class = c("model.specs", 
    "data.frame"), object.info = c(package = "unmarked", func = "distsamp", 
    className = "unmarkedFitDS", fitType = "distsamp"))), unmarkedFitDSO = list(
    distsampOpen = structure(list(model = c("abundance", "recruitment", 
    "recruitment", "growth", "growth", "growth", "survival", 
    "carrying capacity", "detection", "immigration", "scale", 
    "dispersion", "zeroinfl"), item.name = c("lambda", "gamma", 
    "gamma", "gamma", "gamma", "gamma", "omega", "omega", "det", 
    "iota", "scale", "alpha", "psi"), formula.arg = c("lambdaformula", 
    "gammaformula", "gammaformula", "gammaformula", "gammaformula", 
    "gammaformula", "omegaformula", "omegaformula", "pformula", 
    "iotaformula", NA, NA, NA), name = c("Abundance", "Recruitment", 
    "Recruitment", "Growth Rate", "Growth Rate", "Growth Rate", 
    "Apparent Survival", "Carrying Capacity", "Detection", "Immigration", 
    "Hazard-rate(scale)", "Dispersion", "Zero-inflation"), short.name = c("lam", 
    "gamConst", "gamAR", "gamTrend", "gamRicker", "gamGomp", 
    "omega", "omCarCap", "sigma", "iota", "scale", "alpha", "psi"
    ), formula.slot = c("formlist", "formlist", "formlist", "formlist", 
    "formlist", "formlist", "formlist", "formlist", "formlist", 
    "formlist", NA, NA, NA), formlist.name = c("lambdaformula", 
    "gammaformula", "gammaformula", "gammaformula", "gammaformula", 
    "gammaformula", "omegaformula", "omegaformula", "pformula", 
    "iotaformula", NA, NA, NA), formula.slot2 = c("formula1", 
    "formula2", "formula2", "formula2", "formula2", "formula2", 
    "formula3", "formula3", "formula4", "formula5", NA, NA, NA
    ), data = structure(c(1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
    2L, NA, NA, NA), levels = c("site", "site+year"), class = "factor")), row.names = c(NA, 
    -13L), class = c("model.specs", "data.frame"), object.info = c(package = "unmarked", 
    func = "distsampOpen", className = "unmarkedFitDSO", fitType = "distsampOpen"
    )), multmixOpen = structure(list(model = c("recruitment", 
    "recruitment", "growth", "growth", "growth"), item.name = c("gamma", 
    "gamma", "gamma", "gamma", "gamma"), formula.arg = c("gammaformula", 
    "gammaformula", "gammaformula", "gammaformula", "gammaformula"
    ), name = c("Recruitment", "Recruitment", "Growth Rate", 
    "Growth Rate", "Growth Rate"), short.name = c("gamConst", 
    "gamAR", "gamTrend", "gamRicker", "gamGomp"), formula.slot = c("formlist", 
    "formlist", "formlist", "formlist", "formlist"), formlist.name = c("gammaformula", 
    "gammaformula", "gammaformula", "gammaformula", "gammaformula"
    ), formula.slot2 = c("formula2", "formula2", "formula2", 
    "formula2", "formula2"), data = structure(c(1L, 1L, 1L, 1L, 
    1L), levels = "site+year", class = "factor")), row.names = c(NA, 
    -5L), class = c("model.specs", "data.frame"), object.info = c(package = "unmarked", 
    func = "multmixOpen", className = "unmarkedFitDSO", fitType = "multmixOpen"
    ))), unmarkedFitGDR = list(gdistremoval = structure(list(
    model = c("abundance", "dispersion", "availability", "distance", 
    "scale", "removal"), item.name = c("lambda", "alpha", "phi", 
    "dist", "scale", "rem"), formula.arg = c("lambdaformula", 
    NA, "phiformula", "distanceformula", NA, "removalformula"
    ), name = c("Abundance", "Dispersion", "Availability", "Distance", 
    "Hazard-rate (scale)", "Removal"), short.name = c("lambda", 
    "alpha", "phi", "dist", "scale", "rem"), formula.slot = c("formlist", 
    NA, "formlist", "formlist", NA, "formlist"), formlist.name = c("lambdaformula", 
    NA, "phiformula", "distanceformula", NA, "removalformula"
    ), formula.slot2 = c("formula1", NA, "formula2", "formula3", 
    NA, "formula4"), data = structure(c(1L, NA, 2L, 2L, NA, 3L
    ), levels = c("site", "site+year", "site+year+obs"), class = "factor")), row.names = c(NA, 
-6L), class = c("model.specs", "data.frame"), object.info = c(package = "unmarked", 
func = "gdistremoval", className = "unmarkedFitGDR", fitType = "gdistremoval"
))), unmarkedFitGDS = list(gdistsamp = structure(list(model = c("abundance", 
"availability", "detection", "scale", "dispersion", "zeroinfl"
), item.name = c("lambda", "phi", "det", "scale", "alpha", "psi"
), formula.arg = c("lambdaformula", "phiformula", "pformula", 
NA, NA, NA), name = c("Abundance", "Availability", "Detection", 
"Hazard-rate(scale)", "Dispersion", "Zero-inflation"), short.name = c("lambda", 
"phi", "p", "scale", "alpha", "psi"), formula.slot = c("formlist", 
"formlist", "formlist", NA, NA, NA), formlist.name = c("lambdaformula", 
"phiformula", "pformula", NA, NA, NA), formula.slot2 = c("formula1", 
"formula2", "formula3", NA, NA, NA), data = structure(c(1L, 2L, 
3L, NA, NA, NA), levels = c("site", "site+year", "site+year+obs"
), class = "factor")), row.names = c(NA, -6L), class = c("model.specs", 
"data.frame"), object.info = c(package = "unmarked", func = "gdistsamp", 
className = "unmarkedFitGDS", fitType = "gdistsamp"))), unmarkedFitGMM = list(
    gmn = structure(list(model = c("abundance", "availability", 
    "detection", "dispersion", "zeroinfl"), item.name = c("lambda", 
    "phi", "det", "alpha", "psi"), formula.arg = c("lambdaformula", 
    "phiformula", "pformula", NA, NA), name = c("Abundance", 
    "Availability", "Detection", "Dispersion", "Zero-inflation"
    ), short.name = c("lambda", "phi", "p", "alpha", "psi"), 
        formula.slot = c("formlist", "formlist", "formlist", 
        NA, NA), formlist.name = c("lambdaformula", "phiformula", 
        "pformula", NA, NA), formula.slot2 = c("formula", "formula", 
        "formula", NA, NA), data = structure(c(1L, 2L, 3L, NA, 
        NA), levels = c("site", "site+year", "site+year+obs"), class = "factor")), row.names = c(NA, 
    -5L), class = c("model.specs", "data.frame"), object.info = c(package = "unmarked", 
    func = "gmultmix", className = "unmarkedFitGMM", fitType = "gmn"
    ))), unmarkedFitGOccu = list(goccu = structure(list(model = c("occupancy", 
"availability", "detection"), item.name = c("psi", "phi", "det"
), formula.arg = c("psiformula", "phiformula", "pformula"), name = c("Occupancy", 
"Availability", "Detection"), short.name = c("psi", "phi", "p"
), formula.slot = c("formlist", "formlist", "formlist"), formlist.name = c("psiformula", 
"phi", "pformula"), formula.slot2 = c("formula", "formula", "formula"
), data = structure(1:3, levels = c("site", "site+year", "site+year+obs"
), class = "factor")), row.names = c(NA, -3L), class = c("model.specs", 
"data.frame"), object.info = c(package = "unmarked", func = "goccu", 
className = "unmarkedFitGOccu", fitType = "goccu"))), unmarkedFitGPC = list(
    gpcount = structure(list(model = c("abundance", "availability", 
    "detection", "dispersion", "zeroinfl"), item.name = c("lambda", 
    "phi", "det", "alpha", "psi"), formula.arg = c("lambdaformula", 
    "phiformula", "pformula", NA, NA), name = c("Abundance", 
    "Availability", "Detection", "Dispersion", "Zero-inflation"
    ), short.name = c("lambda", "phi", "p", "alpha", "psi"), 
        formula.slot = c("formlist", "formlist", "formlist", 
        NA, NA), formlist.name = c("lambdaformula", "phiformula", 
        "pformula", NA, NA), formula.slot2 = c("formula", "formula", 
        "formula", NA, NA), data = structure(c(1L, 2L, 3L, NA, 
        NA), levels = c("site", "site+year", "site+year+obs"), class = "factor")), row.names = c(NA, 
    -5L), class = c("model.specs", "data.frame"), object.info = c(package = "unmarked", 
    func = "gpcount", className = "unmarkedFitGPC", fitType = "gpcount"
    ))), unmarkedFitMMO = list(multmixOpen = structure(list(model = c("abundance", 
"survival", "carrying capacity", "detection", "immigration", 
"dispersion", "zeroinfl"), item.name = c("lambda", "omega", "omega", 
"det", "iota", "alpha", "psi"), formula.arg = c("lambdaformula", 
"omegaformula", "omegaformula", "pformula", "iotaformula", NA, 
NA), name = c("Abundance", "Apparent Survival", "Carrying Capacity", 
"Detection", "Immigration", "Dispersion", "Zero-inflation"), 
    short.name = c("lam", "omega", "omCarCap", "p", "iota", "alpha", 
    "psi"), formula.slot = c("formlist", "formlist", "formlist", 
    "formlist", "formlist", NA, NA), formlist.name = c("lambdaformula", 
    "omegaformula", "omegaformula", "pformula", "iotaformula", 
    NA, NA), formula.slot2 = c("formula1", "formula3", "formula3", 
    "formula4", "formula5", NA, NA), data = structure(c(1L, 2L, 
    2L, 3L, 2L, NA, NA), levels = c("site", "site+year", "site+year+obs"
    ), class = "factor")), row.names = c(NA, -7L), class = c("model.specs", 
"data.frame"), object.info = c(package = "unmarked", func = "multmixOpen", 
className = "unmarkedFitMMO", fitType = "multmixOpen"))), unmarkedFitMPois = list(
    multinomPois = structure(list(model = c("abundance", "detection"
    ), item.name = c("state", "det"), formula.arg = c("formula2", 
    "formula1"), name = c("Abundance", "Detection"), short.name = c("lambda", 
    "p"), formula.slot = c("formula2", "formula1"), formlist.name = c(NA_character_, 
    NA_character_), formula.slot2 = c(NA_character_, NA_character_
    ), data = structure(1:2, levels = c("site", "site+obs"), class = "factor")), row.names = c(NA, 
    -2L), class = c("model.specs", "data.frame"), object.info = c(package = "unmarked", 
    func = "multinomPois", className = "unmarkedFitMPois", fitType = "multinomPois"
    ))), unmarkedFitNmixTTD = list(nmixTTD = structure(list(model = c("abundance", 
"detection", "dispersion", "shape"), item.name = c("state", "det", 
"alpha", "shape"), formula.arg = c("stateformula", "detformula", 
NA, NA), name = c("Abundance", "Detection", "Dispersion", "Weibull shape"
), short.name = c("lamN", "lamP", "alpha", "k"), formula.slot = c("stateformula", 
"detformula", NA, NA), formlist.name = c(NA_character_, NA_character_, 
NA_character_, NA_character_), formula.slot2 = c("formula2", 
"formula1", NA, NA), data = structure(c(1L, 2L, NA, NA), levels = c("site", 
"site+year+obs"), class = "factor")), row.names = c(NA, -4L), class = c("model.specs", 
"data.frame"), object.info = c(package = "unmarked", func = "nmixTTD", 
className = "unmarkedFitNmixTTD", fitType = "nmixTTD"))), unmarkedFitOccu = list(
    occu = structure(list(model = c("occupancy", "detection"), 
        item.name = c("state", "det"), formula.arg = c("formula2", 
        "formula1"), name = c("Occupancy", "Detection"), short.name = c("psi", 
        "p"), formula.slot = c("formula2", "formula1"), formlist.name = c(NA_character_, 
        NA_character_), formula.slot2 = c(NA_character_, NA_character_
        ), data = structure(1:2, levels = c("site", "site+obs"
        ), class = "factor")), row.names = c(NA, -2L), class = c("model.specs", 
    "data.frame"), object.info = c(package = "unmarked", func = "occu", 
    className = "unmarkedFitOccu", fitType = "occu"))), unmarkedFitOccuCOP = list(
    occuCOP = structure(list(model = c("occupancy", "detection rate"
    ), item.name = c("psi", "lambda"), formula.arg = c("psiformula", 
    "lambdaformula"), name = c("Occupancy probability", "Detection rate"
    ), short.name = c("psi", "lambda"), formula.slot = c("formlist", 
    "formlist"), formlist.name = c("psiformula", "lambdaformula"
    ), formula.slot2 = c("formula", "formula"), data = structure(1:2, levels = c("site", 
    "site+obs"), class = "factor")), row.names = c(NA, -2L), class = c("model.specs", 
    "data.frame"), object.info = c(package = "unmarked", func = "occuCOP", 
    className = "unmarkedFitOccuCOP", fitType = "occuCOP"))), 
    unmarkedFitOccuFP = list(occuFP = structure(list(model = c("occupancy", 
    "detection", "false positive", "certainity"), item.name = c("state", 
    "det", "fp", "b"), formula.arg = c("stateformula", "detformula", 
    "Fpformula", "Bformula"), name = c("Occupancy", "Detection", 
    "false positive", "Pcertain"), short.name = c("psi", "p", 
    "fp", "b"), formula.slot = c("stateformula", "detformula", 
    "FPformula", "Bformula"), formlist.name = c(NA_character_, 
    NA_character_, NA_character_, NA_character_), formula.slot2 = c("formula", 
    "formula", "formula", "formula"), data = structure(c(1L, 
    2L, 2L, 2L), levels = c("site", "site+obs"), class = "factor")), row.names = c(NA, 
    -4L), class = c("model.specs", "data.frame"), object.info = c(package = "unmarked", 
    func = "occuFP", className = "unmarkedFitOccuFP", fitType = "occuFP"
    ))), unmarkedFitOccuMS = list(occuMS = structure(list(model = c("occupancy", 
    "initial occupancy", "detection", "transition"), item.name = c("state", 
    "state", "det", "transition"), formula.arg = c("psiformulas", 
    "psiformulas", "phiformulas", "detformulas"), name = c("Occupancy", 
    "Initial Occupancy", "Detection", "Transition Probabilities"
    ), short.name = c("psi", "psi", "p", "phi"), formula.slot = c("formula", 
    "formula", "formula", "formula"), formlist.name = c(NA_character_, 
    NA_character_, NA_character_, NA_character_), formula.slot2 = c(NA_character_, 
    NA_character_, NA_character_, NA_character_), data = structure(c(1L, 
    1L, 2L, 2L), levels = c("site", "site+year+obs"), class = "factor")), row.names = c(NA, 
    -4L), class = c("model.specs", "data.frame"), object.info = c(package = "unmarked", 
    func = "occuMS", className = "unmarkedFitOccuMS", fitType = "occuMS"
    ))), unmarkedFitOccuMulti = list(occuMulti = structure(list(
        model = c("occupancy", "detection"), item.name = c("state", 
        "det"), formula.arg = c("stateformulas", "detformulas"
        ), name = c("Occupancy", "Detection"), short.name = c("psi", 
        "p"), formula.slot = c("formula", "formula"), formlist.name = c(NA_character_, 
        NA_character_), formula.slot2 = c(NA_character_, NA_character_
        ), data = structure(1:2, levels = c("site", "site+obs"
        ), class = "factor")), row.names = c(NA, -2L), class = c("model.specs", 
    "data.frame"), object.info = c(package = "unmarked", func = "occuMulti", 
    className = "unmarkedFitOccuMulti", fitType = "occuMulti"
    ))), unmarkedFitOccuPEN = list(occu = structure(list(model = c("occupancy", 
    "detection"), item.name = c("state", "det"), formula.arg = c("formula2", 
    "formula1"), name = c("Occupancy", "Detection"), short.name = c("psi", 
    "p"), formula.slot = c("formula2", "formula1"), formlist.name = c(NA_character_, 
    NA_character_), formula.slot2 = c(NA_character_, NA_character_
    ), data = structure(1:2, levels = c("site", "site+obs"), class = "factor")), row.names = c(NA, 
    -2L), class = c("model.specs", "data.frame"), object.info = c(package = "unmarked", 
    func = "occuPEN", className = "unmarkedFitOccuPEN", fitType = "occu"
    ))), unmarkedFitOccuPEN_CV = list(occu = structure(list(model = c("occupancy", 
    "detection"), item.name = c("state", "det"), formula.arg = c("formula2", 
    "formula1"), name = c("Occupancy", "Detection"), short.name = c("psi", 
    "p"), formula.slot = c("formula2", "formula1"), formlist.name = c(NA_character_, 
    NA_character_), formula.slot2 = c(NA_character_, NA_character_
    ), data = structure(1:2, levels = c("site", "site+obs"), class = "factor")), row.names = c(NA, 
    -2L), class = c("model.specs", "data.frame"), object.info = c(package = "unmarked", 
    func = "occuPEN_CV", className = "unmarkedFitOccuPEN_CV", 
    fitType = "occu"))), unmarkedFitOccuRN = list(occuRN = structure(list(
        model = c("abundance", "detection"), item.name = c("state", 
        "det"), formula.arg = c("formula2", "formula1"), name = c("Abundance", 
        "Detection"), short.name = c("lam", "p"), formula.slot = c("formula2", 
        "formula1"), formlist.name = c(NA_character_, NA_character_
        ), formula.slot2 = c(NA_character_, NA_character_), data = structure(1:2, levels = c("site", 
        "site+obs"), class = "factor")), row.names = c(NA, -2L
    ), class = c("model.specs", "data.frame"), object.info = c(package = "unmarked", 
    func = "occuRN", className = "unmarkedFitOccuRN", fitType = "occuRN"
    ))), unmarkedFitOccuTTD = list(occuTTD = structure(list(model = c("occupancy", 
    "detection", "colonization", "extinction", "shape"), item.name = c("psi", 
    "det", "col", "ext", "shape"), formula.arg = c("psiformula", 
    "detformula", "gammaformula", "epsilonformula", NA), name = c("Occupancy", 
    "Detection", "Colonization", "Extinction", "Weibull shape"
    ), short.name = c("psi", "lam", "col", "ext", "k"), formula.slot = c("psiformula", 
    "detformula", "gamformula", "epsformula", NA), formlist.name = c(NA_character_, 
    NA_character_, NA_character_, NA_character_, NA_character_
    ), formula.slot2 = c(NA_character_, NA_character_, NA_character_, 
    NA_character_, NA_character_), data = structure(c(1L, 3L, 
    2L, 2L, NA), levels = c("site", "site+year", "site+year+obs"
    ), class = "factor")), row.names = c(NA, -5L), class = c("model.specs", 
    "data.frame"), object.info = c(package = "unmarked", func = "occuTTD", 
    className = "unmarkedFitOccuTTD", fitType = "occuTTD"))), 
    unmarkedFitPCO = list(pcountOpen = structure(list(model = c("abundance", 
    "recruitment", "recruitment", "growth", "growth", "growth", 
    "survival", "carrying capacity", "detection", "immigration", 
    "dispersion", "zeroinfl"), item.name = c("lambda", "gamma", 
    "gamma", "gamma", "gamma", "gamma", "omega", "omega", "det", 
    "iota", "alpha", "psi"), formula.arg = c("lambdaformula", 
    "gammaformula", "gammaformula", "gammaformula", "gammaformula", 
    "gammaformula", "omegaformula", "omegaformula", "pformula", 
    "iotaformula", NA, NA), name = c("Abundance", "Recruitment", 
    "Recruitment", "Growth Rate", "Growth Rate", "Growth Rate", 
    "Apparent Survival", "Carrying Capacity", "Detection", "Immigration", 
    "Dispersion", "Zero-inflation"), short.name = c("lam", "gamConst", 
    "gamAR", "gamTrend", "gamRicker", "gamGomp", "omega", "omCarCap", 
    "p", "iota", "alpha", "psi"), formula.slot = c("formlist", 
    "formlist", "formlist", "formlist", "formlist", "formlist", 
    "formlist", "formlist", "formlist", "formlist", NA, NA), 
        formlist.name = c("lambdaformula", "gammaformula", "gammaformula", 
        "gammaformula", "gammaformula", "gammaformula", "omegaformula", 
        "omegaformula", "pformula", "iotaformula", NA, NA), formula.slot2 = c("formula1", 
        "formula2", "formula2", "formula2", "formula2", "formula2", 
        "formula3", "formula3", "formula4", "formula5", NA, NA
        ), data = structure(c(1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
        3L, 2L, NA, NA), levels = c("site", "site+year", "site+year+obs"
        ), class = "factor")), row.names = c(NA, -12L), class = c("model.specs", 
    "data.frame"), object.info = c(package = "unmarked", func = "pcountOpen", 
    className = "unmarkedFitPCO", fitType = "pcountOpen"))), 
    unmarkedFitPCount = list(pcount = structure(list(model = c("abundance", 
    "detection", "dispersion", "zeroinfl"), item.name = c("state", 
    "det", "alpha", "psi"), formula.arg = c("formula2", "formula1", 
    NA, NA), name = c("Abundance", "Detection", "Dispersion", 
    "Zero-inflation"), short.name = c("lam", "p", "alpha", "psi"
    ), formula.slot = c("formula2", "formula1", NA, NA), formlist.name = c(NA_character_, 
    NA_character_, NA_character_, NA_character_), formula.slot2 = c(NA_character_, 
    NA_character_, NA_character_, NA_character_), data = structure(c(1L, 
    2L, NA, NA), levels = c("site", "site+obs"), class = "factor")), row.names = c(NA, 
    -4L), class = c("model.specs", "data.frame"), object.info = c(package = "unmarked", 
    func = "pcount", className = "unmarkedFitPCount", fitType = "pcount"
    )))))
