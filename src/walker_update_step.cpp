#include <vector>
#include "cpp11.hpp"
#include "Rmath.h"
#include "local_rng.h"
using namespace cpp11;


void WalkerUpdateWeights(
    integers&, const std::vector<double>&, double, double,
    std::vector<double>&, std::vector<double>&);

void WalkerUpdateClusterPhiTau(
    int, const doubles&, const integers&, double, double, double, double,
    std::vector<double>&, std::vector<double>&);

void WalkerUpdateClusterIdentifiers(
    doubles&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&,
    const std::vector<double>&, std::vector<int>&);

double WalkerUpdateAlpha(
    const std::vector<int>&, double, double, double, int);

double UpdateMuPhi(
    const std::vector<double>&, const std::vector<double>&, double, double, double);

std::vector<double> UpdateCalendarAges(
        int, const doubles&, double, double, const std::vector<int>&, const std::vector<double>&,
        const std::vector<double>&, const doubles&, const doubles&, int, const doubles&, const doubles&);

// Performs one iteration of Walker DP update on all model parameters.
// Updates the following parameters:
// - calendar_ages
// - weights
// - v
// - phi, tau
// - cluster_ids
// - n_clust
// - alpha
// - mu_phi
// Uses the following hyperparameters which are unchanged
// - alpha_shape, alpha_rate, lambda, nu1, nu2, A, B
[[cpp11::register]] list WalkerUpdateStep(
    doubles current_calendar_ages, // current calendar age each for each observation
    doubles current_weight,        // weight per cluster
    doubles current_v,
    integers current_cluster_ids,  // cluster each observation belongs to
    double current_alpha,          // current value of the DPMM concentration parameter
    double current_mu_phi,         // current value of the overall cluster centering
    double alpha_shape,
    double alpha_rate,
    double lambda,
    double nu1,
    double nu2,
    double A,
    double B,
    double w,
    double m,
    doubles c14_determinations,
    doubles c14_sigmas,
    int calcurve_yr_index_offset,
    doubles mucalallyr,
    doubles sigcalallyr) {

  local_rng rng_state;                   // Ensures RNG follows R and R follows after
  int n = current_calendar_ages.size();  // Number of observations
  std::vector<double> u(n);              // Auxiliary variables
  std::vector<double> weight;            // Updated weights
  weight.reserve(2 * current_weight.size());  // Reserve extra space for new clusters
  std::vector<double> v(current_v.begin(), current_v.end());  // Updated v
  std::vector<int> cluster_ids(n);      // Updated cluster_ids
  std::vector<double> phi, tau;         // Updated cluster means and precisions
  std::vector<double> calendar_ages;    // Updated calendar_ages
  double mu_phi, alpha;                 // Updated value of mu_phi and alpha
  double min_u = 1.;                    // Minimum value of auxiliary variables

  using namespace cpp11::literals;
  cpp11::writable::list retlist;

  // Create auxiliary variables u
  for (int k = 0; k < n; k++) {
    u[k] = Rf_runif(0., current_weight[current_cluster_ids[k] - 1]);
    if (u[k] < min_u) min_u = u[k];   // update minimum value of u
  }

  WalkerUpdateWeights(current_cluster_ids, u, min_u, current_alpha, v, weight);

  // Update the cluster means and precisions, introducing new ones for those without observations
  // First set the vectors to the correct size now n_clust has been updated
  phi.resize(weight.size());
  tau.resize(weight.size());
  WalkerUpdateClusterPhiTau(
    weight.size(), current_calendar_ages, current_cluster_ids, current_mu_phi, lambda, nu1, nu2, phi, tau);

  // Now update the cluster id for each observation
  WalkerUpdateClusterIdentifiers(current_calendar_ages, u, weight, phi, tau, cluster_ids);

  alpha = WalkerUpdateAlpha(cluster_ids, current_alpha, alpha_shape, alpha_rate, weight.size());
  mu_phi = UpdateMuPhi(phi, tau, lambda, A, B);

  calendar_ages = UpdateCalendarAges(
    n,
    current_calendar_ages,
    w,
    m,
    cluster_ids,
    phi,
    tau,
    c14_determinations,
    c14_sigmas,
    calcurve_yr_index_offset,
    mucalallyr,
    sigcalallyr);

  // Return the updated parameters
  retlist.push_back({"weight"_nm = weight});
  retlist.push_back({"v"_nm = v});
  retlist.push_back({"cluster_ids"_nm = cluster_ids});
  retlist.push_back({"phi"_nm = phi});
  retlist.push_back({"tau"_nm = tau});
  retlist.push_back({"alpha"_nm = alpha});
  retlist.push_back({"mu_phi"_nm = mu_phi});
  retlist.push_back({"calendar_ages"_nm = calendar_ages});
  return retlist;
}
