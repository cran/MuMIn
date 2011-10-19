\name{Beetle}
\alias{Beetle}
\encoding{utf-8}
\docType{data}
\title{Flour beetle mortality data}
\description{
Mortality of flour beetles (\emph{Tribolium confusu}) due to exposure to gaseous
carbon disulfide CS\eqn{_{2}}{2}, from Bliss (1935)
}

\usage{
data(Beetle)
}

\format{
\code{Beetle} is a data frame with 5 elements.

\describe{
   \item{dose}{The dose of CS\eqn{_{2}}{2} in mg/L}
   \item{n.tested}{Number of beetles tested}
   \item{n.killed}{Number of beetles killed}
   \item{Prop}{A matrix with two columns named \bold{n.killed} and \bold{n.survived}}
   \item{mortality}{Observed mortality rate}

}}

\source{
Bliss C. I. (1935) The calculation of the dosage-mortality curve.
\emph{Annals of Applied Biology}, 22: 134–167
}

\references{
Burnham, K. P. and Anderson, D. R. (2002) \emph{Model selection and multimodel
inference: a practical information-theoretic approach}. 2nd ed.
}

\examples{

# The "Logistic regression example"
# from Burnham & Anderson (2002) chapter 4.11

data(Beetle)

# Fit a global model with all the considered variables
globmod <- glm(Prop ~ dose + I(dose^2) + log(dose) + I(log(dose)^2),
    data=Beetle, family=binomial)


# A logical expression defining the subset of models to use:
# * either log(dose) or dose
# * the quadratic terms can appear only together with linear terms
msubset <- expression(xor(dose, `log(dose)`) & (dose | !`I(dose^2)`)
    & (`log(dose)` | !`I(log(dose)^2)`))


# Table 4.6

# Use 'varying' argument to fit models with different link functions
# Note the use of 'alist' rather than 'list' in order to keep the 'family'
# objects unevaluated
varying.link <- list(family=alist(
    logit = binomial("logit"),
    probit = binomial("probit"),
    cloglog = binomial("cloglog")
    ))

dredge(globmod, subset = msubset, varying = varying.link, rank=AIC)


# Table 4.7 "models justifiable a priori"
(ms3 <- dredge(update(globmod, . ~ dose), fixed = ~dose, rank = AIC,
    varying = varying.link))

mod3 <- get.models(ms3, 1:3)

# Table 4.8. Predicted mortality probability at dose 40.

# a helper functions to calculate confidence intervals on logit scale
logit.ci <- function(p, se, quantile = 2) {
    C. <- exp(quantile * se / (p * (1 - p)))
    p /(p + (1 - p) * c(C., 1/C.))
}

pred <- sapply(mod3, predict, newdata = list(dose = 40), se.fit = TRUE,
    type = "response")
pred <- apply(pred, 1, unlist)[, 1:2] # simplify

# build the table
tab <- rbind(pred, par.avg(pred[,"fit"], pred[,"se.fit"], Weights(ms3),
    revised.var = FALSE)[1:2])
tab <- cbind(
    c(Weights(ms3), NA),
    tab,
    matrix(logit.ci(tab[,"fit"], tab[,"se.fit"], quantile = c(rep(1.96, 3), 2)),
    ncol=2)
    )
colnames(tab) <- c("Akaike weight", "Predicted(40)", "SE", "Lower CI", "Upper CI")
rownames(tab) <- c(as.character(ms3$family), "model averaged")

print(tab, digits = 3, na.print = "")


# Figure 4.3
newdata <- list(dose = seq(min(Beetle$dose), max(Beetle$dose), length.out = 25))
matplot(newdata$dose, sapply(mod3, predict, newdata, type="response"),
    type = "l", xlab = quote(list("Dose of"~ CS[2],(mg/L))),
    ylab = "Mortality"
)



}

\keyword{datasets}
