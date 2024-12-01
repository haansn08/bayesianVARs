# x should be of class bayesianVARs_bvar

as.fsvdraws <- function(x) {
	if (x$sigma_type != "factor") {
		stop("cannot extract fsvdraws object from non-factor model")
	}
	
	ret <- list()
	ret$y <- x$Y
	ret$fac <- aperm(x$fac, c(2,1,3))
	ret$facload <- x$facload
	ret$para <- x$sv_para
	ret$logvar <- x$logvar
	ret$config <- x$config
	class(ret) <- "fsvdraws"
	
	ret
}

irf <- function(x, shock, ahead=8) {
	M <- ncol(x$Y)
	ret <- irf_cpp(
		x$PHI,
		shock,
		ahead
	)
	dimnames(ret)[[2]] <- dimnames(x$Y)[[2]]
	ret
}
