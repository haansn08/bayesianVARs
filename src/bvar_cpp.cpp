#include <RcppArmadillo.h>
#include <stochvol.h>
#include <progress.hpp>
#include <Rcpp/Benchmark/Timer.h>
#include "sample_coefficients.h"

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
List bvar_cpp(const arma::mat Y,
              const arma::mat X,
              const int M,
              const int T,
              const int K,
              const int draws,
              const int burnin,
              const int intercept,
              const arma::vec priorIntercept,
              arma::mat PHI,
              arma::mat PHI0, //??? no const
              const List priorPHI_in,
              const List priorL_in,
              arma::mat L,
              const bool SV,
              const arma::vec priorHomoscedastic,
              const List sv_spec,
              arma::mat h,
              arma::mat sv_para,
              const arma::imat i_mat,
              const arma::ivec i_vec,
              const bool progressbar){

//-------------------------Preliminaries--------------------------------------//

  const double n = K*M; // number of VAR coefficients // ??? without intercept
  arma::mat PHI_diff(K+intercept,M); // will hold PHI - PHI0
  //const arma::uvec i_ol = arma::find(i_vec > 0); // indicator for ownlags
  //const arma::uvec i_cl = arma::find(i_vec < 0); // indicator for crosslags
  // int n_ol = i_ol.size(); // nr. of ownlags
  //const int n_cl = i_cl.size(); // nr. of crosslags
  const arma::uvec i_ocl= arma::find(i_vec != 0.); // indicator for all coefficients except intercept
  const arma::uvec i_i= arma::find(i_vec == 0.); // indicator for intercepts

  //////////////////
  // indicators for coefficients without intercept
  const arma::ivec i_vec_small=i_vec(i_ocl);
  const arma::uvec i_ol = arma::find(i_vec_small > 0);
  const arma::uvec i_cl = arma::find(i_vec_small < 0);
  const int n_ol = i_ol.size(); // nr. of ownlags
  const int n_cl = i_cl.size(); // nr. of crosslags
  /////////////////

//--------------------Initialization of hyperparameters-----------------------//

//---- PHI
  std::string priorPHI = priorPHI_in["prior"];
  // V_i holds prior variances (without intercepts)
  arma::vec V_i(n);

  arma::vec V_i_long(n+M*intercept); // ??? V_i plus intercept prior variances
  V_i_long(i_i) = priorIntercept; // in case of HM prior, these will be scaled later (probably not???)

  if(priorPHI == "normal"){
    arma::vec V_i_in = priorPHI_in["V_i"];
    //in case of 'normal' V_i is fixed at user specified values
    V_i = V_i_in;
  }

  //---- DL prior on PHI

//  int n_groups;
  int n_groups = priorPHI_in["n_groups"];
//  IntegerVector groups_in;
  IntegerVector groups_in = priorPHI_in["groups"];
//  if(priorPHI == "R2D2" || priorPHI== "DL"){
    //n_groups = priorPHI_in["n_groups"];
    //groups_in = priorPHI_in["groups"];
//  }
  arma::ivec groups(groups_in.begin(), groups_in.length(), false);

  arma::vec psi(n);psi.fill(1.0);
  arma::vec zeta(n_groups); zeta.fill(10);
  arma::vec theta(n);theta.fill(1/static_cast<double>(n));

  NumericVector DL_a_in = priorPHI_in["DL_a"];
  arma::vec DL_a(DL_a_in.begin(), DL_a_in.length(), false);
//  arma::vec DL_a(n_groups);
//  bool DL_hyper;
  bool DL_hyper = priorPHI_in["DL_hyper"];
//  NumericVector a_vec_in;
  NumericVector a_vec_in = priorPHI_in["a_vec"];
//  NumericMatrix prep2_in(n_groups, 1000);
  NumericMatrix prep2_in = priorPHI_in["prep2"];
//  NumericVector prep1_in;
  NumericVector prep1_in = priorPHI_in["prep1"];
  if(priorPHI == "DL" ){
//    DL_hyper = priorPHI_in["DL_hyper"];
    V_i = psi % theta % theta * zeta(0) * zeta(0);
//    arma::vec DL_a_in = priorPHI_in["DL_a"];
//    DL_a = DL_a_in;
//    if(DL_hyper == true){
//      a_vec_in = priorPHI_in["a_vec"];
//      prep2_in = wrap(priorPHI_in["prep2"]);
//      prep1_in = priorPHI_in["prep1"];
//      }
  }
  arma::vec a_vec(a_vec_in.begin(), a_vec_in.length(), false);
  arma::mat prep2(prep2_in.begin(), prep2_in.nrow(), prep2_in.ncol(), false);
  arma::vec prep1(prep1_in.begin(), prep1_in.length(), false);

  //----R2D2 on PHI
  vec xi(n_groups, fill::ones);
  arma::vec theta_r2d2(n); theta_r2d2.fill(1/static_cast<double>(n));
  vec zeta_r2d2(n_groups); zeta_r2d2.fill(10);
  arma::vec psi_r2d2(n); psi_r2d2.fill(1/static_cast<double>(n));

//  bool R2D2_hyper;
  bool R2D2_hyper = priorPHI_in["R2D2_hyper"];
//  NumericVector api_in;
  NumericVector api_in = priorPHI_in["R2D2_api"];
//  NumericVector b_r2d2_in;
  NumericVector b_r2d2_in = priorPHI_in["R2D2_b"];
//  NumericVector api_vec_in;
  NumericVector api_vec_in = priorPHI_in["api_vec"];
//  NumericVector b_r2d2_vec_in;
  NumericVector b_r2d2_vec_in = priorPHI_in["b_vec"];
  if(priorPHI == "R2D2"){

//    api_in = priorPHI_in["R2D2_api"];
//    b_r2d2_in = priorPHI_in["R2D2_b"];

//    R2D2_hyper = priorPHI_in["R2D2_hyper"];
//    if(R2D2_hyper==true){
//          b_r2d2_vec_in = priorPHI_in["b_vec"];
//          api_vec_in = priorPHI_in["api_vec"];
//    }

    V_i = psi_r2d2%theta_r2d2*zeta_r2d2(0)/2;

  }
  arma::vec api(api_in.begin(), api_in.length(), false);
  arma::vec b_r2d2(b_r2d2_in.begin(), b_r2d2_in.length(), false);
  arma::vec b_r2d2_vec(b_r2d2_vec_in.begin(), b_r2d2_vec_in.length(), false);
  arma::vec api_vec(api_vec_in.begin(), api_vec_in.length(), false);

  //---- SSVS on PHI
  arma::vec tau_0;
  arma::vec tau_1;
  double SSVS_s_a;
  double SSVS_s_b;
  if(priorPHI == "SSVS"){
    arma::vec tau_0_in = priorPHI_in["SSVS_tau0"];
    tau_0 = tau_0_in;
    arma::vec tau_1_in = priorPHI_in["SSVS_tau1"];
    tau_1 = tau_1_in;
    double SSVS_s_a_in = priorPHI_in["SSVS_s_a"];
    SSVS_s_a = SSVS_s_a_in;
    double SSVS_s_b_in = priorPHI_in["SSVS_s_b"];
    SSVS_s_b = SSVS_s_b_in;

    V_i = tau_0 % tau_0;
  }
  arma::vec gammas(n, fill::zeros);
  arma::vec p_i(n); p_i.fill(0.5);

  //---- HMP on PHI
  double lambda_1=0.04;
  double lambda_2=0.0016;
  NumericVector V_i_prep_in;
  NumericVector s_r_1_in;
  NumericVector s_r_2_in;
  if(priorPHI == "HMP"){
     V_i_prep_in = priorPHI_in["V_i_prep"];
     s_r_1_in = priorPHI_in["lambda_1"];
     s_r_2_in = priorPHI_in["lambda_2"];
     }
  arma::vec V_i_prep(V_i_prep_in.begin(), V_i_prep_in.length(), false);
  if(priorPHI == "HMP"){
    //V_i_long(i_i) = priorIntercept % V_i_prep(i_i);
    arma::vec V_i_prep_small = V_i_prep(i_ocl);
    V_i(i_ol) = lambda_1*V_i_prep_small(i_ol);
    V_i(i_cl) = lambda_2*V_i_prep_small(i_cl);
  }
  arma::vec s_r_1(s_r_1_in.begin(), s_r_1_in.length(), false);
  arma::vec s_r_2(s_r_2_in.begin(), s_r_2_in.length(), false);

  //Fill V_i_long with remaining prior variances
  V_i_long(i_ocl) = V_i; // ???
  arma::mat V_prior = arma::reshape(V_i_long, K + intercept, M);

//-------------- L
  std::string priorL = priorL_in["prior"];

  uvec L_upper_indices = trimatu_ind( size(L),  1);
  arma::vec l = L(L_upper_indices);
  const double n_L = l.size();
  arma::vec V_i_L(n_L);

  if(priorL == "normal"){
    arma::vec V_i_L_in = priorL_in["V_i"];
    V_i_L = V_i_L_in;
  }
  //---- DL prior on L
  arma::vec psi_L(n_L); psi_L.fill(1.0);
  double zeta_L=10;
  arma::vec theta_L(n_L); theta_L.fill(1/static_cast<double>(n_L));

  bool DL_L_hyper;
  double DL_b;
  NumericVector b_vec_in;
  NumericVector prep2_L_in;
  NumericVector prep1_L_in;
  if(priorL == "DL"){
    double DL_b_in = priorL_in["DL_b"];
    DL_b = DL_b_in;

    DL_L_hyper = priorL_in["DL_hyper"];
    if(DL_L_hyper){
      b_vec_in = priorL_in["b_vec"];
      prep2_L_in = priorL_in["prep2"];
      prep1_L_in = priorL_in["prep1"];
    }

    V_i_L= psi_L% theta_L % theta_L * zeta_L * zeta_L;
  }
  arma::vec b_vec(b_vec_in.begin(), b_vec_in.length(), false);
  arma::vec prep2_L(prep2_L_in.begin(), prep2_L_in.length(), false);
  arma::vec prep1_L(prep1_L_in.begin(), prep1_L_in.length(), false);

  //----R2D2 on L

  double xi_L = 1;
  arma::vec theta_L_r2d2(n_L); theta_L_r2d2.fill(1/static_cast<double>(n_L));
  double zeta_L_r2d2 = 10;
  arma::vec psi_L_r2d2(n_L); psi_L_r2d2.fill(1/static_cast<double>(n_L));

  bool R2D2_L_hyper;
  NumericVector api_vec_L_in;
  NumericVector b_vec_L_r2d2_in;
  double b_L_r2d2;
  double api_L;
  if(priorL == "R2D2"){
    double b_L_r2d2_in = priorL_in["R2D2_b"];
    b_L_r2d2 = b_L_r2d2_in;
    double api_L_in = priorL_in["R2D2_api"];
    api_L = api_L_in;

    R2D2_L_hyper = priorL_in["R2D2_hyper"];
    if(R2D2_L_hyper){
      api_vec_L_in = priorL_in["api_vec"];
      b_vec_L_r2d2_in = priorL_in["b_vec"];
    }
    V_i_L = psi_L_r2d2 % theta_L_r2d2 * zeta_L_r2d2 /2.;
  }
  arma::vec api_vec_L(api_vec_L_in.begin(), api_vec_L_in.length(),false);
  arma::vec b_vec_L_r2d2(b_vec_L_r2d2_in.begin(), b_vec_L_r2d2_in.length(), false);

  //---- SSVS on L
  arma::vec tau_0_L;
  arma::vec tau_1_L;
  double SSVS_s_a_L;
  double SSVS_s_b_L;
  if(priorL == "SSVS"){
    arma::vec tau_0_L_in = priorL_in["SSVS_tau0"];
    tau_0_L = tau_0_L_in;
    arma::vec tau_1_L_in = priorL_in["SSVS_tau1"];
    tau_1_L = tau_1_L_in;
    double SSVS_s_a_L_in = priorL_in["SSVS_s_a"];
    SSVS_s_a_L = SSVS_s_a_L_in;
    double SSVS_s_b_L_in = priorL_in["SSVS_s_b"];
    SSVS_s_b_L = SSVS_s_b_L_in;

    V_i_L = tau_0_L % tau_0_L;
  }
  arma::vec gammas_L(n_L, fill::zeros);
  arma::vec p_i_L(n_L); p_i_L.fill(0.5);

  //---- HMP on L
  double lambda_3 = 0.001;
  NumericVector s_r_3_in;
  if(priorL == "HMP"){
    s_r_3_in = priorL_in["lambda_3"];
    arma::vec V_i_L_tmp(n_L); V_i_L_tmp.fill(1.0);
    V_i_L= lambda_3*V_i_L_tmp;
  }
  arma::vec s_r_3(s_r_3_in.begin(), s_r_3_in.length(), false);

  //-----------------------------SV-settings----------------------------------//
  // Import sv_spec
  NumericVector sv_priormu = sv_spec["priormu"] ;
  NumericVector sv_priorphi= sv_spec["priorphi"];
  NumericVector sv_priorsigma2 = sv_spec["priorsigma2"] ;
  NumericVector sv_offset = sv_spec["sv_offset"];

  // create prior specification object for the update_*_sv functions
  using stochvol::PriorSpec;
  const PriorSpec prior_spec = {
    PriorSpec::Latent0(),  // stationary prior distribution on h0
    PriorSpec::Mu(PriorSpec::Normal(sv_priormu[0], sv_priormu[1])),  // N(mu,sd)
    PriorSpec::Phi(PriorSpec::Beta(sv_priorphi[0], sv_priorphi[1])),  //
    PriorSpec::Sigma2(PriorSpec::Gamma(sv_priorsigma2[0], sv_priorsigma2[1]))  //
  };  // stochvol would even offen more :heavy-tailed, leverage, regression

  //expert settings: these are the same settings as the default settings for the
  //idiosyncratic variances from package factorstochvol
  //(https://github.com/gregorkastner/factorstochvol, accessed 2021-11-12)
  using stochvol::ExpertSpec_FastSV;
  const ExpertSpec_FastSV expert_sv {
    true,  // interweave
    stochvol::Parameterization::CENTERED,  // centered baseline always
    1e-8,  // B011inv,
    1e-12,  //B022inv,
    2,  // MHsteps,
    ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE,  // independece proposal for sigma
    -1,  // MHcontrol unused for independence prior for sigma
    ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL  // immediately reject (mu,phi,sigma) if proposed phi is outside (-1, 1)
  };

  // initialize mixing indicators
  arma::umat mixind(T,M); mixind.fill(5);
  //initialize d_sqrt
  arma::mat d_sqrt = arma::exp(h/2);


  //------------------------------------STORAGE-------------------------------//

  arma::cube sv_latent_draws(draws, T, M);
  arma::cube sv_para_draws(draws, 4,M);
  arma::cube PHI_draws(draws, (K+intercept), M); // ??? K + intercept
  arma::cube L_draws(draws, M, M);

  int phi_hyperparameter_size(0);
  if(priorPHI == "DL" ){
    phi_hyperparameter_size += 2*n_groups + 2*n; // a + zeta + n(theta + psi)
  }else if(priorPHI == "R2D2"){
    phi_hyperparameter_size += 4*n_groups + 2*n; // (b+)xi + zeta + n(theta + psi)
  }else if(priorPHI == "SSVS"){
    phi_hyperparameter_size += 2*n; // n(gammas + p_i)
  }else if(priorPHI == "HMP"){
    phi_hyperparameter_size += 2; // lambda_1 + lambda_2
  }
  arma::mat phi_hyperparameter_draws(draws, phi_hyperparameter_size);

  int l_hyperparameter_size(0);
  if(priorL == "DL" ){
    l_hyperparameter_size += 2. + 2*n_L;
  }else if(priorL == "R2D2"){
    l_hyperparameter_size += 4. + 2*n_L; // b + api +xi + zeta + n(theta + psi)
  }else if(priorL == "SSVS"){
    l_hyperparameter_size += 2*n_L;
  }else if(priorL == "HMP"){
    l_hyperparameter_size += 1;
  }
  arma::mat l_hyperparameter_draws(draws, l_hyperparameter_size);

  //indicator vector needed for DL and R2D2 on L
  arma::uvec l_ind(n_L);
  for(int i=0; i<n_L; ++i){
    l_ind(i)=i;
  }
  //-----------------------------------SAMPLER--------------------------------//

  const int tot = draws + burnin;
  // Initialize progressbar
  Progress p(tot, progressbar);
  Timer timer;
  timer.step("start");
  for(int rep = 0; rep < tot; rep++){

    //----1) Draw PHI (reduced form VAR coefficients)

    try{
      sample_PHI(PHI, PHI0, Y, X, L, d_sqrt, V_prior, K, M, false);
    } catch(...){
      ::Rf_error("Couldn't sample PHI in rep %i.", rep);
    }

    if(priorPHI == "DL" || priorPHI == "R2D2"){
      PHI_diff = PHI - PHI0;
      // if regularization gets extreme, often there appear zeros (numerical issue)
      // coefficients must not be zero, otherwise problems with do_rgig1
      // anyhow, a realized value of a continous pd cannot be exactly zero
      for (int ii = 0; ii<K; ii++) {
        for (int jj = 0; jj<M; jj++){
          if(PHI_diff(ii,jj) == 0) {

            if(R::rbinom( 1, 0.5 )==0){
              PHI(ii,jj) = PHI0(ii,jj) + 1e-100;//eps(ii,jj);//1e-100;
              PHI_diff(ii,jj) = 1e-100;
            }else{
              PHI(ii,jj) = PHI0(ii,jj) - 1e-100;//eps(ii,jj);//1e-100;
              PHI_diff(ii,jj) = -1e-100;
            }
          }else if(PHI_diff(ii,jj) < 1e-100 && PHI_diff(ii,jj) > 0){ //eps(ii,jj)
            PHI(ii,jj) = PHI0(ii,jj) + 1e-100;//eps(ii,jj);
            PHI_diff(ii,jj) = 1e-100;
          }else if (PHI_diff(ii,jj) > -1e-100 && PHI_diff(ii,jj) < 0){ //-eps(ii,jj)
            PHI(ii,jj) = PHI0(ii,jj) - 1e-100;//eps(ii,jj);
            PHI_diff(ii,jj) = -1e-100;
          }
        }
      }
    }else{
      PHI_diff = PHI - PHI0;
    }

    //----2) Sample hyperparameters of hierarchical priors (prior variances V_i)

    if(priorPHI == "DL" || priorPHI== "R2D2"){

      arma::ivec::iterator g;
      int j=0;

      for(g=groups.begin(); g!=groups.end(); ++g){
        arma::uvec ind = arma::find(i_vec_small == *g);
        arma::uvec indplus = arma::find(i_vec == *g);

        if(priorPHI=="DL"){
          sample_V_i_DL(V_i, PHI_diff(indplus), DL_a(j) , a_vec, prep1,
                        prep2.row(j), zeta(j), psi, theta, ind, DL_hyper);
        }else if(priorPHI=="R2D2"){
          sample_V_i_R2D2(V_i, PHI_diff(indplus), api(j), api_vec, zeta_r2d2(j),
                          psi_r2d2, theta_r2d2, xi(j), b_r2d2(j), b_r2d2_vec, ind,
                          R2D2_hyper);
        }


        j += 1;
      }

    }else if(priorPHI == "SSVS"){

      if(rep > 0.1*burnin){
        sample_V_i_SSVS(V_i, gammas, p_i, PHI_diff(i_ocl), tau_0, tau_1,
                        SSVS_s_a, SSVS_s_b);

      }

    }else if(priorPHI == "HMP"){

      if(rep > 0.1*burnin){
        sample_V_i_HMP(lambda_1, lambda_2, V_i, s_r_1(0), s_r_1(1), s_r_2(0),
                       s_r_2(1), PHI_diff(i_ocl), V_i_prep(i_ocl), n_ol, n_cl,
                       i_ol, i_cl);
      }
    }

    V_i_long(i_ocl) = V_i;

    V_prior = reshape(V_i_long, K+intercept, M);

    //----3) Draw Sigma_t = inv(L)'*D_t*inv(L), where L is upper unitriangular
    //       and D diagonal
    //----3a) Draw free off-diagonal elements in L

    arma::mat resid = Y - X*PHI;

    try{
      sample_L(L, resid, V_i_L, d_sqrt);
    }
    catch(...){
      ::Rf_error("Couldn't sample L in rep %i.", rep);
    }

    if(priorL == "DL" || priorL == "R2D2"){

      for (int jj = 0; jj<M; jj++){
        for (int ii = 0; ii<jj; ii++) {
          if(L(ii,jj) == 0) {
            if(R::rbinom( 1, 0.5 )==0){
              L(ii,jj) = 1e-100;
            }else{
              L(ii,jj) = -1e-100;
            }
          }else
            if(L(ii,jj) < 1e-100 && L(ii,jj) > 0){
              L(ii,jj) = 1e-100;
            }else if (L(ii,jj) > -1e-100 && L(ii,jj) < 0){
              L(ii,jj) = -1e-100;
            }
        }
      }
    }

    //----3b) Draw hyperparameters of hierarchical priors
    l = L(L_upper_indices);
    if(priorL == "DL" ){


      try{
        sample_V_i_DL(V_i_L, l, DL_b , b_vec, prep1_L, prep2_L, zeta_L, psi_L,
                      theta_L, l_ind, DL_L_hyper);
      } catch (...) {
        ::Rf_error("Couldn't sample V_i_L (DL prior)  in run %i", rep);

      }

    }else if(priorL == "R2D2"){

      sample_V_i_R2D2(V_i_L, l, api_L, api_vec_L, zeta_L_r2d2, psi_L_r2d2,
                      theta_L_r2d2, xi_L, b_L_r2d2, b_vec_L_r2d2, l_ind, R2D2_L_hyper );

    }else if(priorL == "SSVS"){

      sample_V_i_SSVS(V_i_L, gammas_L, p_i_L, l, tau_0_L, tau_1_L, SSVS_s_a_L, SSVS_s_b_L);

    }else if(priorL == "HMP"){

      sample_V_i_L_HMP(lambda_3, V_i_L, s_r_3(0), s_r_3(1), l);
    }

    //----3c) Draw elements of D_t
    //        in case of SV use package stochvol
    arma::mat str_resid = resid*L; // structural (orthogonalized) residuals

    if(SV == false ){
      for(int i =0; i<M; i++){
        double s_p = priorHomoscedastic(i,1) + 0.5*accu(square(str_resid.col(i)));
        double d_i = 1. / R::rgamma(priorHomoscedastic(i,0)+T/2, 1./s_p);
        d_sqrt.col(i).fill(sqrt(d_i));
        h.col(i).fill(log(d_i));
      }
    }else if(SV == true){
      //const arma::mat resid_norm = log(square(str_resid) + sv_offset); // + 1e-40offset??
      for(int j=0; j < M; j++){
        const arma::mat resid_norm = log(square(str_resid.col(j)) + sv_offset(j));
        arma::vec h_j  = h.unsafe_col(j);  // unsafe_col reuses memory, h will automatically be overwritten
        arma::uvec mixind_j = mixind.unsafe_col(j);
        double mu = sv_para(0,j),
          phi = sv_para(1,j),
          sigma = sv_para(2,j),
          h0_j = sv_para(3,j);
        stochvol::update_fast_sv(resid_norm, mu, phi, sigma, h0_j, h_j, mixind_j, prior_spec, expert_sv); //resid_norm.col(j)
        sv_para.col(j) = arma::colvec({mu, phi, sigma, h0_j});
      }
      d_sqrt = exp(h/2);
    }

    //-------Store draws after burnin
    if(rep >= burnin){

      PHI_draws.row(rep-burnin) = PHI;
      L_draws.row(rep-burnin) = L;
      sv_latent_draws.row(rep-burnin) = h;
      sv_para_draws.row((rep-burnin)) = sv_para;

      if(priorPHI == "DL" ){

        phi_hyperparameter_draws(rep-burnin, span(0,(n_groups-1))) = zeta;
        phi_hyperparameter_draws(rep-burnin, span(n_groups,(n_groups+n-1))) = trans(psi.as_col());
        phi_hyperparameter_draws(rep-burnin, span((n_groups+n),(n_groups+2*n-1)));
        phi_hyperparameter_draws(rep-burnin, span((n_groups+2*n),phi_hyperparameter_size-1.)) = DL_a;

      }else if(priorPHI == "R2D2"){

        phi_hyperparameter_draws(rep-burnin, span(0,(n_groups-1))) = zeta_r2d2 ;
        phi_hyperparameter_draws(rep-burnin, span(n_groups,(n_groups+n-1))) = trans(psi_r2d2.as_col());
        phi_hyperparameter_draws(rep-burnin, span((n_groups+n),(n_groups+2*n-1))) = trans(theta_r2d2.as_col());
        phi_hyperparameter_draws(rep-burnin, span((n_groups+2*n),phi_hyperparameter_size-1.-2*n_groups)) = xi ;
        phi_hyperparameter_draws(rep-burnin, span((2*n_groups+2*n ),phi_hyperparameter_size -1.-n_groups)) = b_r2d2;
        phi_hyperparameter_draws(rep-burnin, span((3*n_groups+2*n ),phi_hyperparameter_size-1.)) = api;

      }else if(priorPHI == "SSVS"){

        phi_hyperparameter_draws(rep-burnin, span(0, (n-1.))) = gammas;
        phi_hyperparameter_draws(rep-burnin, span(n, (phi_hyperparameter_size-1.))) = p_i;

      }else if(priorPHI == "HMP"){
        phi_hyperparameter_draws(rep-burnin, 0) = lambda_1;
        phi_hyperparameter_draws(rep-burnin, 1) = lambda_2;
      }

      if(priorL == "DL"){

        l_hyperparameter_draws(rep-burnin, 0) = zeta_L;
        l_hyperparameter_draws(rep-burnin, span(1,n_L)) = trans(psi_L.as_col());
        l_hyperparameter_draws(rep-burnin, span(n_L+1.,2*n_L)) = trans(theta_L.as_col());
        l_hyperparameter_draws(rep-burnin, l_hyperparameter_size-1.) = DL_b;

      }else if(priorL == "R2D2"){

        l_hyperparameter_draws(rep-burnin, 0) = zeta_L_r2d2 ;
        l_hyperparameter_draws(rep-burnin, span(1,(n_L))) = trans(psi_L_r2d2.as_col());
        l_hyperparameter_draws(rep-burnin, span(n_L+1.,(2*n_L))) = trans(theta_L_r2d2.as_col());
        l_hyperparameter_draws(rep-burnin, l_hyperparameter_size-3.) = xi_L ;
        l_hyperparameter_draws(rep-burnin, l_hyperparameter_size-2.) = b_L_r2d2;
        l_hyperparameter_draws(rep-burnin, l_hyperparameter_size-1.) = api_L;

      }else if(priorL == "SSVS"){

        l_hyperparameter_draws(rep-burnin, span(0, (n_L-1.))) = gammas_L;
        l_hyperparameter_draws(rep-burnin, span(n_L, (l_hyperparameter_size-1.))) = p_i_L;

      }else if(priorL == "HMP"){

        l_hyperparameter_draws(rep-burnin, 0) = lambda_3;
      }
    }

    p.increment();
  }
  timer.step("end");
  NumericVector time(timer);

  List out = List::create(
    Named("PHI") = PHI_draws,
    Named("L") = L_draws,
    Named("sv_latent") = sv_latent_draws,
    Named("sv_para") = sv_para_draws,
    Named("phi_hyperparameter") = phi_hyperparameter_draws,
    Named("l_hyperparameter") = l_hyperparameter_draws,
    Named("bench") = time,
    //Named("i_ocl") = i_ocl, // ???
    //Named("i_i") = i_i, // ???
    //Named("s_r_3") = s_r_3, // ???
    //Named("tau0") = tau_0, // ???
    //Named("tau1") = tau_1,// ???
    Named("V_i_long") = V_i_long, // ???
    Named("V_i_L") = V_i_L,
    //Named("a_vec") = a_vec, // ???
    Named("groups") = groups,
    Named("api") = api,
    Named("DL_a") = DL_a,
    Named("prep1") = prep1,
    Named("prep2") = prep2,// ???
    Named("PHI0")=PHI0,
    Named("i_vec_small") = i_vec_small,
    Named("PHI_diff") = PHI_diff,
    Named("n_groups") = n_groups,
    Named("groups") = groups

  );

  return out;
}
