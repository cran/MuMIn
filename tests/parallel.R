if(MuMIn:::.parallelPkgCheck(quiet = TRUE)) {
	clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
	clust <- try(makeCluster(getOption("cl.cores", 2), type = clusterType))
	if(inherits(clust, "cluster")) {
	
		library(MuMIn)
		library(nlme)

		data(Orthodont, package = "nlme")
		#Orthodont <- Orthodont[sample.int(nrow(Orthodont), size = 64,
			#replace = TRUE), ]
		Orthodont$rand1 <- runif(nrow(Orthodont))
		Orthodont$rand2 <- runif(nrow(Orthodont))
		clusterExport(clust, "Orthodont")
		clusterCall(clust, "library", "nlme", character.only = TRUE)		

		# fm2 <- lmer(log(distance) ~ rand*Sex*age + (1|Subject) + (1|Sex),
		#	data = Orthodont, REML=FALSE)
		fm2 <- lme(log(distance) ~ rand1*Sex*age + rand2, ~ 1|Subject / Sex,
			data = Orthodont, method = "ML")
		print(system.time(pdd1 <- dredge(fm2, cluster = FALSE)))
		print(system.time(pddc <- dredge(fm2, cluster = clust)))
		print(system.time(dd1 <- dredge(fm2)))

		print(pddc)
		print(pdd1)
		print(dd1)
		
		#print(all.equal(pddc, dd1))
		ma1 <- model.avg(pdd1, beta = "none")
		ma0 <- model.avg(pddc)

		if(!isTRUE(test <- all.equal(ma1$avg.model, ma0$avg.model))) {
			print(test)
			warning("'ma1' and 'ma0' are not equal")
		}
		if(!isTRUE(test <- all.equal(ma1$summary, ma0$summary))) {
			print(test)
			warning("'ma1' and 'ma0' are not equal")
		}
		
        if(!(identical(c(pddc), c(pdd1)) && identical(c(pdd1), c(dd1)))) {
		    warning("results of 'dredge' and 'pdredge' are not equal")
			print(all.equal(c(pddc), c(pdd1)))
			print(all.equal(c(pdd1), c(dd1)))
		}

		stopCluster(clust)

	# suppressPackageStartupMessages(library(spdep))
	# suppressMessages(example(NY_data, echo = FALSE))
	# esar1f <- spautolm(Z ~ PEXPOSURE * PCTAGE65P + PCTOWNHOME,
	# data=nydata, listw=listw_NY, family="SAR", method="full", verbose=FALSE)
	# clusterCall(clust, "library", "spdep", character.only = TRUE)
	# clusterExport(clust, "listw_NY", "nydata")
	# options(warn=1)

	# varying <- list(family = list("CAR", "SAR"), method=list("Matrix_J", "full"))

	# dd <- dredge(esar1f, m.lim=c(0, 1),  fixed=~PEXPOSURE, varying = varying, trace=FALSE)

	} else # if(inherits(clust, "try-error"))
	message("Could not set up the cluster")

}

#system.time(pdredge(fm2, cluster = clust))
#system.time(pdredge(fm2, cluster = F))
#system.time(dredge(fm2))
