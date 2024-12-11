// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bvar_cpp
List bvar_cpp(const arma::mat& Y, const arma::mat& X, const int& M, const int& T, const int& K, const int& draws, const int& burnin, const int& thin, const std::string& tvp_keep, const int& intercept, const arma::vec priorIntercept, arma::mat& PHI0, const List priorPHI_in, const List priorSigma_in, const List Rstartvals_in, const arma::imat& i_mat, const arma::ivec& i_vec, const bool& progressbar, const double& PHI_tol, const double& L_tol, const bool& huge);
RcppExport SEXP _bayesianVARs_bvar_cpp(SEXP YSEXP, SEXP XSEXP, SEXP MSEXP, SEXP TSEXP, SEXP KSEXP, SEXP drawsSEXP, SEXP burninSEXP, SEXP thinSEXP, SEXP tvp_keepSEXP, SEXP interceptSEXP, SEXP priorInterceptSEXP, SEXP PHI0SEXP, SEXP priorPHI_inSEXP, SEXP priorSigma_inSEXP, SEXP Rstartvals_inSEXP, SEXP i_matSEXP, SEXP i_vecSEXP, SEXP progressbarSEXP, SEXP PHI_tolSEXP, SEXP L_tolSEXP, SEXP hugeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int& >::type T(TSEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type draws(drawsSEXP);
    Rcpp::traits::input_parameter< const int& >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< const int& >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type tvp_keep(tvp_keepSEXP);
    Rcpp::traits::input_parameter< const int& >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type priorIntercept(priorInterceptSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type PHI0(PHI0SEXP);
    Rcpp::traits::input_parameter< const List >::type priorPHI_in(priorPHI_inSEXP);
    Rcpp::traits::input_parameter< const List >::type priorSigma_in(priorSigma_inSEXP);
    Rcpp::traits::input_parameter< const List >::type Rstartvals_in(Rstartvals_inSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type i_mat(i_matSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type i_vec(i_vecSEXP);
    Rcpp::traits::input_parameter< const bool& >::type progressbar(progressbarSEXP);
    Rcpp::traits::input_parameter< const double& >::type PHI_tol(PHI_tolSEXP);
    Rcpp::traits::input_parameter< const double& >::type L_tol(L_tolSEXP);
    Rcpp::traits::input_parameter< const bool& >::type huge(hugeSEXP);
    rcpp_result_gen = Rcpp::wrap(bvar_cpp(Y, X, M, T, K, draws, burnin, thin, tvp_keep, intercept, priorIntercept, PHI0, priorPHI_in, priorSigma_in, Rstartvals_in, i_mat, i_vec, progressbar, PHI_tol, L_tol, huge));
    return rcpp_result_gen;
END_RCPP
}
// my_gig
NumericMatrix my_gig(int n, NumericVector lambda, NumericVector chi, NumericVector psi);
RcppExport SEXP _bayesianVARs_my_gig(SEXP nSEXP, SEXP lambdaSEXP, SEXP chiSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type chi(chiSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(my_gig(n, lambda, chi, psi));
    return rcpp_result_gen;
END_RCPP
}
// obtain_restrictable_matrices
Rcpp::List obtain_restrictable_matrices(const arma::cube& reduced_coefficients, const arma::cube& factor_loadings, const arma::mat& U_vecs, const arma::mat& logvar_T, const bool include_B0, const bool include_structural_coeff, const bool include_long_run_ir);
RcppExport SEXP _bayesianVARs_obtain_restrictable_matrices(SEXP reduced_coefficientsSEXP, SEXP factor_loadingsSEXP, SEXP U_vecsSEXP, SEXP logvar_TSEXP, SEXP include_B0SEXP, SEXP include_structural_coeffSEXP, SEXP include_long_run_irSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type reduced_coefficients(reduced_coefficientsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type factor_loadings(factor_loadingsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type U_vecs(U_vecsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type logvar_T(logvar_TSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_B0(include_B0SEXP);
    Rcpp::traits::input_parameter< const bool >::type include_structural_coeff(include_structural_coeffSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_long_run_ir(include_long_run_irSEXP);
    rcpp_result_gen = Rcpp::wrap(obtain_restrictable_matrices(reduced_coefficients, factor_loadings, U_vecs, logvar_T, include_B0, include_structural_coeff, include_long_run_ir));
    return rcpp_result_gen;
END_RCPP
}
// irf_cpp
arma::cube irf_cpp(const arma::cube& coefficients, const arma::cube& factor_loadings, const arma::mat& U_vecs, const arma::colvec& shock, const arma::uword ahead);
RcppExport SEXP _bayesianVARs_irf_cpp(SEXP coefficientsSEXP, SEXP factor_loadingsSEXP, SEXP U_vecsSEXP, SEXP shockSEXP, SEXP aheadSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type coefficients(coefficientsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type factor_loadings(factor_loadingsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type U_vecs(U_vecsSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type shock(shockSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type ahead(aheadSEXP);
    rcpp_result_gen = Rcpp::wrap(irf_cpp(coefficients, factor_loadings, U_vecs, shock, ahead));
    return rcpp_result_gen;
END_RCPP
}
// sample_PHI_cholesky
arma::mat sample_PHI_cholesky(const arma::mat PHI, const arma::mat& PHI_prior, const arma::mat& Y, const arma::mat& X, const arma::mat& U, const arma::mat& d_sqrt, const arma::mat& V_prior);
RcppExport SEXP _bayesianVARs_sample_PHI_cholesky(SEXP PHISEXP, SEXP PHI_priorSEXP, SEXP YSEXP, SEXP XSEXP, SEXP USEXP, SEXP d_sqrtSEXP, SEXP V_priorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type PHI(PHISEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type PHI_prior(PHI_priorSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type d_sqrt(d_sqrtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type V_prior(V_priorSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_PHI_cholesky(PHI, PHI_prior, Y, X, U, d_sqrt, V_prior));
    return rcpp_result_gen;
END_RCPP
}
// dmvnrm_arma_fast
arma::vec dmvnrm_arma_fast(arma::mat const& x, arma::rowvec const& mean, arma::mat const& sigma, bool const logd);
RcppExport SEXP _bayesianVARs_dmvnrm_arma_fast(SEXP xSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP logdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat const& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::rowvec const& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::mat const& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool const >::type logd(logdSEXP);
    rcpp_result_gen = Rcpp::wrap(dmvnrm_arma_fast(x, mean, sigma, logd));
    return rcpp_result_gen;
END_RCPP
}
// predh
arma::cube predh(const arma::mat& logvar_T, const arma::ivec& ahead, const int& each, const arma::mat& sv_mu, const arma::mat& sv_phi, const arma::mat& sv_sigma);
RcppExport SEXP _bayesianVARs_predh(SEXP logvar_TSEXP, SEXP aheadSEXP, SEXP eachSEXP, SEXP sv_muSEXP, SEXP sv_phiSEXP, SEXP sv_sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type logvar_T(logvar_TSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ahead(aheadSEXP);
    Rcpp::traits::input_parameter< const int& >::type each(eachSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sv_mu(sv_muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sv_phi(sv_phiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sv_sigma(sv_sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(predh(logvar_T, ahead, each, sv_mu, sv_phi, sv_sigma));
    return rcpp_result_gen;
END_RCPP
}
// out_of_sample
Rcpp::List out_of_sample(const int& each, const arma::rowvec& X_T_plus_1, const arma::cube& PHI, const arma::mat& U, const arma::cube& facload, const arma::mat& logvar_T, const arma::ivec& ahead, const arma::mat& sv_mu, const arma::mat& sv_phi, const arma::mat& sv_sigma, const arma::uvec& sv_indicator, const bool& factor, const bool& LPL, const arma::mat& Y_obs, const bool& LPL_subset, const arma::urowvec& VoI, const bool& simulate_predictive);
RcppExport SEXP _bayesianVARs_out_of_sample(SEXP eachSEXP, SEXP X_T_plus_1SEXP, SEXP PHISEXP, SEXP USEXP, SEXP facloadSEXP, SEXP logvar_TSEXP, SEXP aheadSEXP, SEXP sv_muSEXP, SEXP sv_phiSEXP, SEXP sv_sigmaSEXP, SEXP sv_indicatorSEXP, SEXP factorSEXP, SEXP LPLSEXP, SEXP Y_obsSEXP, SEXP LPL_subsetSEXP, SEXP VoISEXP, SEXP simulate_predictiveSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type each(eachSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type X_T_plus_1(X_T_plus_1SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type PHI(PHISEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type facload(facloadSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type logvar_T(logvar_TSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type ahead(aheadSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sv_mu(sv_muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sv_phi(sv_phiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sv_sigma(sv_sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type sv_indicator(sv_indicatorSEXP);
    Rcpp::traits::input_parameter< const bool& >::type factor(factorSEXP);
    Rcpp::traits::input_parameter< const bool& >::type LPL(LPLSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Y_obs(Y_obsSEXP);
    Rcpp::traits::input_parameter< const bool& >::type LPL_subset(LPL_subsetSEXP);
    Rcpp::traits::input_parameter< const arma::urowvec& >::type VoI(VoISEXP);
    Rcpp::traits::input_parameter< const bool& >::type simulate_predictive(simulate_predictiveSEXP);
    rcpp_result_gen = Rcpp::wrap(out_of_sample(each, X_T_plus_1, PHI, U, facload, logvar_T, ahead, sv_mu, sv_phi, sv_sigma, sv_indicator, factor, LPL, Y_obs, LPL_subset, VoI, simulate_predictive));
    return rcpp_result_gen;
END_RCPP
}
// insample
arma::cube insample(const arma::mat& X, const arma::cube& PHI, const arma::mat& U, const arma::cube& facload, const arma::cube& logvar, const bool& prediction, const bool& factor);
RcppExport SEXP _bayesianVARs_insample(SEXP XSEXP, SEXP PHISEXP, SEXP USEXP, SEXP facloadSEXP, SEXP logvarSEXP, SEXP predictionSEXP, SEXP factorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type PHI(PHISEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type facload(facloadSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type logvar(logvarSEXP);
    Rcpp::traits::input_parameter< const bool& >::type prediction(predictionSEXP);
    Rcpp::traits::input_parameter< const bool& >::type factor(factorSEXP);
    rcpp_result_gen = Rcpp::wrap(insample(X, PHI, U, facload, logvar, prediction, factor));
    return rcpp_result_gen;
END_RCPP
}
// vcov_cpp
arma::cube vcov_cpp(const bool& factor, const arma::cube& facload, const arma::cube& logvar, const arma::mat& U, const int& M, const int& factors);
RcppExport SEXP _bayesianVARs_vcov_cpp(SEXP factorSEXP, SEXP facloadSEXP, SEXP logvarSEXP, SEXP USEXP, SEXP MSEXP, SEXP factorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const bool& >::type factor(factorSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type facload(facloadSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type logvar(logvarSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type U(USEXP);
    Rcpp::traits::input_parameter< const int& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const int& >::type factors(factorsSEXP);
    rcpp_result_gen = Rcpp::wrap(vcov_cpp(factor, facload, logvar, U, M, factors));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayesianVARs_bvar_cpp", (DL_FUNC) &_bayesianVARs_bvar_cpp, 21},
    {"_bayesianVARs_my_gig", (DL_FUNC) &_bayesianVARs_my_gig, 4},
    {"_bayesianVARs_obtain_restrictable_matrices", (DL_FUNC) &_bayesianVARs_obtain_restrictable_matrices, 7},
    {"_bayesianVARs_irf_cpp", (DL_FUNC) &_bayesianVARs_irf_cpp, 5},
    {"_bayesianVARs_sample_PHI_cholesky", (DL_FUNC) &_bayesianVARs_sample_PHI_cholesky, 7},
    {"_bayesianVARs_dmvnrm_arma_fast", (DL_FUNC) &_bayesianVARs_dmvnrm_arma_fast, 4},
    {"_bayesianVARs_predh", (DL_FUNC) &_bayesianVARs_predh, 6},
    {"_bayesianVARs_out_of_sample", (DL_FUNC) &_bayesianVARs_out_of_sample, 17},
    {"_bayesianVARs_insample", (DL_FUNC) &_bayesianVARs_insample, 7},
    {"_bayesianVARs_vcov_cpp", (DL_FUNC) &_bayesianVARs_vcov_cpp, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayesianVARs(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
