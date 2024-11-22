> class(mod)
[1] "bayesianVARs_bvar"

# construct a svdraws object because we don't have one :(
object <- list()
object$svlm <- FALSE
object$meanmodel <- "none"
object$thinning <- list(para = 1, latent = 1, time = "all")

object$para <- coda::mcmc.list()
object$latent <- coda::mcmc.list()
number_of_sv_processes <- dim(mod$sv_para)[2]
#for the "factor" setting
#number_of_sv_processes = M + number of factors
for (i in seq_len(number_of_sv_processes)) {
	sv_para_draws <- cbind(t(mod$sv_para[,i,]), Inf, 0)
	colnames(sv_para_draws) <- c("mu", "phi", "sigma", "nu", "rho")
	object$para[[i]] <- coda::as.mcmc(sv_para_draws)
	
	sv_latent <- t(mod$logvar[,i,])
	object$latent[[i]] <- coda::as.mcmc(sv_latent)
}

class(object) <- c("svdraws")
