#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::cube shock_propagating_predict_cpp(
	const arma::cube& coefficients, //rows: lagged variables + intercept, columns: variables, slices: draws
	const arma::rowvec predictors, //predictors at t=T, columns: lagged variables + 1
	const arma::uword ahead, //how far to predict ahead
	const arma::uword each, //how many times per posterior sample we draw from the predictive distr.
	const arma::cube& facload, //rows: variables, cols: factors, slices: draws
	const arma::cube& predicted_idih, //rows: variables, cols: draws * each, slices: ahead
	const arma::cube& predicted_factorh //rows: factors, cols: draws * each, slices: ahead
) {
	const uword n_variables = coefficients.n_cols;
	const uword n_posterior_draws = coefficients.n_slices;
	arma::cube predictions(ahead, n_variables, each * n_posterior_draws, arma::fill::none);
	
	// trace out each*n_posterior_draws paths
	arma::cube current_predictors(each, predictors.n_cols, n_posterior_draws);
	// all paths start with the given predictors
	for (uword r = 0; r < n_posterior_draws; r++) {
		current_predictors.slice(r).each_row() = predictors;
	}
	for (uword t = 0; t < ahead; t++) {
		arma::mat idiosyncratic_errors(predicted_idih.n_rows, predicted_idih.n_cols, arma::fill::none);
		idiosyncratic_errors.imbue(R::norm_rand);
		idiosyncratic_errors %= arma::exp(predicted_idih.slice(t)/2);
		
		arma::mat factor_draws(predicted_factorh.n_rows, predicted_factorh.n_cols, arma::fill::none);
		factor_draws.imbue(R::norm_rand);
		factor_draws %= arma::exp(predicted_factorh.slice(t)/2);
		
		for(uword r = 0; r < n_posterior_draws; r++) {
			// predict new values for each path
			// rows: each, cols: predicted new variables
			arma::mat new_predictions = current_predictors.slice(r) * coefficients.slice(r);
			
			// add random shocks
			new_predictions += (facload.slice(r) * factor_draws.cols(r, r+each-1)).t();
			new_predictions += idiosyncratic_errors.cols(r, r+each-1).t();
			
			for (uword j = 0; j < each; j++) {
				// store the predictions in a safe place
				predictions.slice(j + r * each).row(t) = new_predictions.row(j);
			}
			
			// make predictions the new predictors keeping the intercept
			current_predictors.slice(r).head_cols(n_variables) = new_predictions;
		}		
	}
	
	return predictions;
}
