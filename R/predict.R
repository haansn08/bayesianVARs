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

shock_propagating_predict <- function(x, ahead=1, each=1) {
	M <- ncol(x$Y)
	predictors <- as.numeric(x$datamat[nrow(x$datamat), 1:(x$lags*M)])
  	if(x$intercept) predictors <- c(predictors, 1)
  	
  	#predict factor specific and idiosyncratic variances for the factor model
  	h_pred <- factorstochvol::predh(as.fsvdraws(x), ahead=seq_len(ahead), each=each)
	
	ret <- shock_propagating_predict_cpp(
		x$PHI,
		predictors,
		ahead,
		each,
		x$facload,
		h_pred$idih,
		h_pred$factorh
	)
	dimnames(ret)[[2]] <- dimnames(x$Y)[[2]]
	ret
}
