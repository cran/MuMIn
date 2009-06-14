`coeffs.default` <-
function(model) { model$coefficients}

`coeffs.glmer` <-
function(model) { ret <- model@fixef; names(ret) <- model@cnames$.fixed; ret}

`coeffs.gls` <-
function (model) return(summary(model)$coefficients)

`coeffs.lme` <-
function(model) { model$coefficients$fixed}

`coeffs.lmer` <-
function(model) { ret <- model@fixef; names(ret) <- model@cnames$.fixed; ret}

`coeffs.mer` <-
function(model) { return(model@fixef)}

`coeffs` <-
function (model) UseMethod("coeffs")

`coeffs.spautolm` <-
function(model) { model$fit$coefficients}
