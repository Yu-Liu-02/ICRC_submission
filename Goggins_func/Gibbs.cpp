// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

inline int Z_at(double t_e, int k_i, int m1, const NumericVector &tgrid) {
  if (k_i <= m1)
    return (t_e >= tgrid[k_i - 1]) ? 1 : 0;
  return 0;
}

// [[Rcpp::export]]
List gibbs_sample_k_with_w_cpp(NumericVector alpha, double beta,
                               NumericVector w, NumericVector time,
                               IntegerVector status, NumericMatrix X,
                               IntegerVector k_init, NumericVector tgrid,
                               List feasible_list, int burn, int thin, int keep,
                               int seed) {

  const int n = time.size();
  const int p = X.ncol();
  const int m1 = tgrid.size();
  const int mK = w.size();

  IntegerVector k = clone(k_init);

  std::vector<int> ev;
  ev.reserve(n);
  for (int i = 0; i < n; ++i)
    if (status[i] == 1)
      ev.push_back(i);
  std::sort(ev.begin(), ev.end(),
            [&](int a, int b) { return time[a] < time[b]; });
  const int E = (int)ev.size();

  std::vector<double> tE(E);
  for (int e = 0; e < E; ++e)
    tE[e] = time[ev[e]];

  // PRECOMPUTE X_i^T * alpha for all patients
  std::vector<double> eta_base(n, 0.0);
  for (int j = 0; j < n; ++j) {
    double tmp = 0.0;
    for (int d = 0; d < p; ++d) {
      tmp += alpha[d] * X(j, d);
    }
    eta_base[j] = tmp;
  }

  std::vector<double> D(E, 0.0);
  for (int e = 0; e < E; ++e) {
    double te = tE[e];
    double sum = 0.0;
    for (int j = 0; j < n; ++j) {
      if (time[j] >= te)
        sum += std::exp(eta_base[j] + beta * Z_at(te, k[j], m1, tgrid));
    }
    D[e] = sum;
  }

  IntegerMatrix draws(keep, n);
  const int total_iters = burn + keep * thin;
  int saved = 0;

  for (int iter = 1; iter <= total_iters; ++iter) {

    if (iter % 50 == 1 && iter > 1) {
      for (int e = 0; e < E; ++e) {
        double te = tE[e];
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
          if (time[j] >= te)
            sum += std::exp(eta_base[j] + beta * Z_at(te, k[j], m1, tgrid));
        }
        D[e] = sum;
      }
    }

    for (int i = 0; i < n; ++i) {
      IntegerVector cand = feasible_list[i];
      const int C = cand.size();
      if (C <= 1)
        continue;

      const int k_old = k[i];
      int eMax = -1;
      for (int e = 0; e < E; ++e) {
        if (tE[e] <= time[i])
          eMax = e;
        else
          break;
      }

      std::vector<double> lps(C, 0.0);
      std::vector<double> exp_old;
      if (eMax >= 0)
        exp_old.resize(eMax + 1);

      if (eMax < 0) {
        for (int cc = 0; cc < C; ++cc) {
          int kc = cand[cc];
          lps[cc] = std::log(std::max(w[kc - 1], 1e-300));
        }
      } else {
        double E0 = std::exp(eta_base[i]);
        double E1 = std::exp(eta_base[i] + beta);

        std::vector<double> V0(eMax + 1, 0.0);
        std::vector<double> V1(eMax + 1, 0.0);

        for (int e = 0; e <= eMax; ++e) {
          int zold = Z_at(tE[e], k_old, m1, tgrid);
          exp_old[e] = (zold == 1) ? E1 : E0;

          double D_base = std::max(D[e] - exp_old[e], 1e-300);
          double D_new0 = std::max(D_base + E0, 1e-300);
          double D_new1 = std::max(D_base + E1, 1e-300);

          V0[e] = std::log(D[e]) - std::log(D_new0);
          V1[e] = std::log(D[e]) - std::log(D_new1);

          if (ev[e] == i) {
            double eta_old = eta_base[i] + beta * zold;
            V0[e] += eta_base[i] - eta_old;
            V1[e] += (eta_base[i] + beta) - eta_old;
          }
        }

        std::vector<double> sum_V0(eMax + 2, 0.0);
        std::vector<double> sum_V1(eMax + 2, 0.0);
        for (int e = 0; e <= eMax; ++e) {
          sum_V0[e + 1] = sum_V0[e] + V0[e];
          sum_V1[e + 1] = sum_V1[e] + V1[e];
        }

        for (int cc = 0; cc < C; ++cc) {
          int kc = cand[cc];
          double lp = std::log(std::max(w[kc - 1], 1e-300));

          double flip_time = (kc <= m1) ? tgrid[kc - 1] : INFINITY;
          auto it =
              std::lower_bound(tE.begin(), tE.begin() + eMax + 1, flip_time);
          int e_split = std::distance(tE.begin(), it);

          lp += sum_V0[e_split] + (sum_V1[eMax + 1] - sum_V1[e_split]);
          lps[cc] = lp;
        }
      }

      double mx = lps[0];
      for (int cc = 1; cc < C; ++cc)
        if (lps[cc] > mx)
          mx = lps[cc];

      std::vector<double> p_prob(C);
      double s = 0.0;
      for (int cc = 0; cc < C; ++cc) {
        p_prob[cc] = std::exp(lps[cc] - mx);
        s += p_prob[cc];
      }

      double u = R::runif(0.0, 1.0) * s;
      double cdf = 0.0;
      int pick = C - 1;
      for (int cc = 0; cc < C; ++cc) {
        cdf += p_prob[cc];
        if (u <= cdf) {
          pick = cc;
          break;
        }
      }

      int k_new = cand[pick];

      if (k_new != k_old && eMax >= 0) {
        for (int e = 0; e <= eMax; ++e) {
          int znew = Z_at(tE[e], k_new, m1, tgrid);
          D[e] = std::max(
              D[e] - exp_old[e] + std::exp(eta_base[i] + beta * znew), 1e-300);
        }
        k[i] = k_new;
      }
    }

    if (iter > burn && ((iter - burn) % thin == 0)) {
      for (int j = 0; j < n; ++j)
        draws(saved, j) = k[j];
      saved++;
      if (saved >= keep)
        break;
    }
  }

  return List::create(_["draws"] = draws, _["k_last"] = k);
}

// [[Rcpp::export]]
double mean_ploglik_ab_cpp(NumericVector alpha, double beta, NumericVector time,
                           IntegerVector status, NumericMatrix X,
                           IntegerMatrix draws, NumericVector tgrid) {

  const int n = time.size();
  const int p = X.ncol();
  const int m1 = tgrid.size();
  const int B = draws.nrow();

  std::vector<int> ev;
  ev.reserve(n);
  for (int i = 0; i < n; ++i)
    if (status[i] == 1)
      ev.push_back(i);
  if (ev.size() == 0)
    return 0.0;
  std::sort(ev.begin(), ev.end(),
            [&](int a, int b) { return time[a] < time[b]; });
  const int E = ev.size();

  std::vector<std::vector<int>> Rset(E);
  for (int e = 0; e < E; ++e) {
    double t = time[ev[e]];
    for (int j = 0; j < n; ++j) {
      if (time[j] >= t)
        Rset[e].push_back(j);
    }
  }

  // PRECOMPUTE base risk for Z=0 and Z=1 across all covariates
  std::vector<double> eta0(n), eta1(n);
  for (int j = 0; j < n; ++j) {
    double tmp = 0.0;
    for (int d = 0; d < p; ++d)
      tmp += alpha[d] * X(j, d);
    eta0[j] = tmp;
    eta1[j] = tmp + beta;
  }

  double total_ll = 0.0;

  for (int b = 0; b < B; ++b) {
    IntegerVector k = draws(b, _);
    double ll_b = 0.0;

    for (int e = 0; e < E; ++e) {
      double t = time[ev[e]];
      int i_ev = ev[e];

      int Zi = (k[i_ev] <= m1 && t >= tgrid[k[i_ev] - 1]) ? 1 : 0;
      double eta_i = (Zi == 1) ? eta1[i_ev] : eta0[i_ev];

      double max_eta = -INFINITY;
      for (int j : Rset[e]) {
        double val = (k[j] <= m1 && t >= tgrid[k[j] - 1]) ? eta1[j] : eta0[j];
        if (val > max_eta)
          max_eta = val;
      }

      double sum_exp = 0.0;
      for (int j : Rset[e]) {
        double val = (k[j] <= m1 && t >= tgrid[k[j] - 1]) ? eta1[j] : eta0[j];
        sum_exp += std::exp(val - max_eta);
      }

      ll_b += eta_i - (max_eta + std::log(sum_exp));
    }
    total_ll += ll_b;
  }
  return total_ll / B;
}