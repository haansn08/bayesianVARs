#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::cube shock_propagating_predict(
	const arma::cube& coefficients_posterior_draws, //rows: lagged variables + intercept, columns: variables, slices: draws
	const arma::rowvec predictors, //predictors at t=T, columns: lagged variables + 1
	const arma::uword ahead, //how far to predict ahead
	const arma::uword each //how many times per posterior sample we draw from the predictive distr.
) {
	const uword n_variables = coefficients_posterior_draws.n_cols;
	const uword n_posterior_draws = coefficients_posterior_draws.n_slices;
	arma::cube predictions(ahead, n_variables, each * n_posterior_draws, arma::fill::none);
	
	// trace out each*n_posterior_draws paths
	arma::cube current_predictors(each, predictors.n_cols, coefficients_posterior_draws.n_slices);
	// all paths start with the given predictors
	for (uword r = 0; r < n_posterior_draws; r++) {
		current_predictors.slice(r).each_row() = predictors;
	}
	for (uword t = 0; t < ahead; t++) {
		for(uword r = 0; r < n_posterior_draws; r++) {
			// predict new values for each path
			// rows: each, cols: predicted new variables
			arma::mat new_predictions = current_predictors.slice(r) * coefficients_posterior_draws.slice(r);
			
			// add random shocks
			// TODO: let the shocks follow the SV process
			arma::mat shocks(new_predictions.n_rows, new_predictions.n_cols, arma::fill::none);
			shocks.imbue(R::norm_rand);
			new_predictions += shocks;
			
			// store the predictions in a safe place
			for (uword j = 0; j < each; j++) {
				predictions.slice(j + r * each).row(t) = new_predictions.row(j);
			}
			
			// make predictions the new predictors keeping the intercept
			current_predictors.slice(r).head_cols(n_variables) = new_predictions;
		}		
	}
	
	return predictions;
}
