if (MuMIn:::testStart("survival")) {

    bladder1 <- bladder[bladder$enum < 5, ]

    fmcph <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) + cluster(id), bladder1)

    r.squared.coxph <- function(object, ...) {
        logtest <- -2 * (object$loglik[1L] - object$loglik[2L])
        c(rsq = 1 - exp(-logtest/object$n), maxrsq = 1 - exp(2 * object$loglik[1L]/object$n))
    }

    getAllTerms(fmcph)
    coef(fmcph)

    ms <- dredge(fmcph, fixed=c("cluster(id)", "strata(enum)"), 
        extra = list(R2 = "r.squared.coxph"))

    summary(model.avg(ms))

    # BUG in survival
    if(! "logLik.coxph.null" %in% methods("logLik"))
        registerS3method("logLik", "coxph.null", survival:::logLik.coxph.null)

    summary(model.avg(ms[1:10]))

    fits <- get.models(ms, delta < 5)
    summary(model.avg(fits))

    ####

    lung <- na.omit(lung)
    fm <- coxph(Surv(time, status) ~ ph.ecog + tt(age), data=lung,
         tt=function(x,t,...) pspline(x + t/365.25))
    ma <- model.avg(dredge(fm))
    coef(ma)
    coefTable(ma)
    ####

    fmsrvrg <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='weibull',
        scale = 1)

    summary(model.avg(dredge(fmsrvrg), delta  < 4))

    fmsrvrg2 <- survreg(Surv(futime, fustat) ~ ecog.ps + rx, ovarian, dist='weibull')

    fmsrvrg3 <- survreg(Surv(time, status) ~ ph.ecog + age + strata(sex), lung,
          na.action = "na.omit")

    coefTable(fmsrvrg)
    coefTable(fmsrvrg2)
    coefTable(fmsrvrg3)
}