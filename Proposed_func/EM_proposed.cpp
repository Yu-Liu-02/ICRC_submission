// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// g_type: 0 = indicator 1(s >= t),  1 = ReLU max(s-t, 0)
inline double g_cpp(double s, double t, int g_type) {
  if (g_type == 0)
    return (s >= t) ? 1.0 : 0.0;
  return std::max(s - t, 0.0);
}

// Helper: X %*% v for n x p matrix X and length-p vector v -> length-n output
static void matvec(const NumericMatrix &X, const NumericVector &v,
                   std::vector<double> &out) {
  int n = X.nrow(), p = X.ncol();
  std::fill(out.begin(), out.end(), 0.0);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < p; j++)
      out[i] += X(i, j) * v[j];
}

// ============================================================
// 1. calculate_p_cpp
// Returns n x (m1+1) matrix
// ============================================================
// [[Rcpp::export]]
NumericMatrix calculate_p_cpp(
    NumericMatrix X,     // n x p
    NumericVector Y,     // n  (observed times)
    IntegerVector delta, // n  (event indicator)
    NumericVector R, // n  (right endpoint of IC interval; Inf = right-censored)
    NumericVector gamma,     // p  (IC covariate effects)
    NumericVector alpha,     // p  (RC covariate effects)
    double beta,             // scalar (single intermediate-event effect)
    NumericVector cumlambda, // m1+1 (cumulative IC baseline hazard)
    NumericVector h,         // m2   (RC baseline hazard increments)
    NumericVector t,         // m1   (IC time grid, sorted)
    NumericVector s,         // m2   (RC event times, sorted)
    IntegerVector k_start,   // n    (1-indexed, IC interval start)
    IntegerVector k_end,     // n    (1-indexed, IC interval end)
    int g_type = 0           // 0 = indicator, 1 = ReLU
) {
  int n = X.nrow(), p = X.ncol();
  int m1 = t.size(), m2 = s.size();

  NumericMatrix P(n, m1 + 1);
  std::fill(P.begin(), P.end(), 0.0);

  double exp_beta = std::exp(beta);

  std::vector<double> xg(n, 0.0), xa(n, 0.0);
  matvec(X, gamma, xg);
  matvec(X, alpha, xa);
  std::vector<double> exp_gammaX(n), exp_alphaX(n);
  for (int i = 0; i < n; i++) {
    exp_gammaX[i] = std::exp(xg[i]);
    exp_alphaX[i] = std::exp(xa[i]);
  }

  // Prefix sums of h: h_cum[l+1] = sum_{j=0}^{l} h[j]
  std::vector<double> h_cum(m2 + 1, 0.0);
  for (int l = 0; l < m2; l++)
    h_cum[l + 1] = h_cum[l] + h[l];

  for (int i = 0; i < n; i++) {
    double eGX = exp_gammaX[i];
    double eAX = exp_alphaX[i];

    // n_risk: number of s values <= Y[i]
    int n_risk = (int)(std::upper_bound(s.begin(), s.end(), Y[i]) - s.begin());

    // spec_idx: where s[l] == Y[i] (only needed for events)
    int spec_idx = -1;
    if (delta[i] == 1) {
      for (int l = 0; l < m2; l++)
        if (s[l] == Y[i]) {
          spec_idx = l;
          break;
        }
    }

    // --- IC intervals ---
    if (k_start[i] <= k_end[i]) {
      for (int k = k_start[i] - 1; k <= k_end[i] - 1; k++) {
        double tk = t[k];
        double p_IC =
            std::exp(-cumlambda[k] * eGX) - std::exp(-cumlambda[k + 1] * eGX);

        // cut: first l where s[l] >= tk  (g transitions from 0 to nonzero here)
        int cut = (int)(std::lower_bound(s.begin(), s.begin() + n_risk, tk) -
                        s.begin());

        double sum_h;
        if (g_type == 0) {
          // Indicator: O(1) prefix-sum split
          sum_h = h_cum[cut] + exp_beta * (h_cum[n_risk] - h_cum[cut]);
        } else {
          // ReLU: g=0 side uses prefix sum; g>0 side loops (s[l]-tk > 0)
          sum_h = h_cum[cut];
          for (int l = cut; l < n_risk; l++)
            sum_h += h[l] * std::exp(beta * (s[l] - tk));
        }
        double temp_k = std::exp(-eAX * sum_h);

        double p_RC_k = temp_k;
        if (delta[i] == 1 && spec_idx >= 0) {
          double g_Y_tk = (g_type == 0) ? ((Y[i] >= tk) ? 1.0 : 0.0)
                                        : std::max(Y[i] - tk, 0.0);
          p_RC_k *= h[spec_idx] * eAX * std::exp(beta * g_Y_tk);
        }

        P(i, k) = p_IC * p_RC_k;
      }
    }

    // --- RC / right-censored column (index m1): t_k unknown so g = 0 ---
    if (k_start[i] > k_end[i] || std::isinf(R[i])) {
      double p_IC = std::exp(-cumlambda[m1] * eGX);
      double sum_h = h_cum[n_risk]; // g=0: just prefix sum
      double temp = std::exp(-eAX * sum_h);

      double p_RC = temp;
      if (delta[i] == 1 && spec_idx >= 0)
        p_RC *= h[spec_idx] * eAX;

      P(i, m1) = p_IC * p_RC;
    }
  }

  return P;
}

// ============================================================
// 2. calculate_w_cpp  (no g — unchanged)
// Returns n x m1 weight matrix
// ============================================================
// [[Rcpp::export]]
NumericMatrix calculate_w_cpp(NumericMatrix X,       // n x p
                              NumericVector gamma,   // p
                              NumericVector lambda,  // m1
                              NumericMatrix q,       // n x (m1+1)
                              IntegerVector k_start, // n (1-indexed)
                              IntegerVector k_end    // n (1-indexed)
) {
  int n = X.nrow(), p = X.ncol();
  int m1 = lambda.size();

  NumericMatrix w(n, m1);
  std::fill(w.begin(), w.end(), 0.0);

  std::vector<double> xg(n, 0.0);
  matvec(X, gamma, xg);
  std::vector<double> exp_gammaX(n);
  for (int i = 0; i < n; i++)
    exp_gammaX[i] = std::exp(xg[i]);

  for (int i = 0; i < n; i++) {
    if (k_start[i] > k_end[i])
      continue;

    double eGX = exp_gammaX[i];
    double cumq = 0.0;

    for (int k = k_start[i] - 1; k <= k_end[i] - 1; k++) {
      double lk = lambda[k];
      double el = std::exp(-lk * eGX);

      if (el < 1.0)
        w(i, k) = q(i, k) * lk * eGX / (1.0 - el);

      if (k > k_start[i] - 1)
        w(i, k) += cumq * lk * eGX;

      cumq += q(i, k);
    }
  }

  return w;
}

// ============================================================
// 3. update_h_cpp
// Returns m2 vector h
// ============================================================
// [[Rcpp::export]]
NumericVector update_h_cpp(NumericMatrix X,     // n x p
                           NumericVector Y,     // n
                           IntegerVector delta, // n
                           NumericVector s,     // m2 (sorted)
                           NumericVector t,     // m1 (sorted)
                           NumericMatrix q,     // n x (m1+1)
                           NumericVector alpha, // p
                           double beta,         // scalar
                           int g_type = 0       // 0 = indicator, 1 = ReLU
) {
  int n = X.nrow(), p = X.ncol();
  int m1 = t.size(), m2 = s.size();

  NumericVector h(m2, 0.0);

  double exp_beta = std::exp(beta);

  std::vector<double> xa(n, 0.0);
  matvec(X, alpha, xa);
  std::vector<double> exp_alphaX(n);
  for (int i = 0; i < n; i++)
    exp_alphaX[i] = std::exp(xa[i]);

  // Row-wise prefix sums of q (IC columns 0..m1-1)
  std::vector<double> cumq_ic((size_t)n * (m1 + 1), 0.0);
  for (int i = 0; i < n; i++)
    for (int k = 0; k < m1; k++)
      cumq_ic[i * (m1 + 1) + k + 1] = cumq_ic[i * (m1 + 1) + k] + q(i, k);

  for (int l = 0; l < m2; l++) {
    double sl = s[l];

    double numer = 0.0;
    for (int i = 0; i < n; i++)
      if (delta[i] == 1 && Y[i] == sl)
        numer += 1.0;
    if (numer == 0.0)
      continue;

    // cl: number of t[k] <= sl  (g > 0 for these k under both g types)
    int cl = (int)(std::upper_bound(t.begin(), t.end(), sl) - t.begin());

    // Pre-compute exp(beta * g(sl, t[k])) only for k < cl (g > 0 side)
    std::vector<double> exp_bg(cl);
    for (int k = 0; k < cl; k++) {
      double gv = (g_type == 0) ? 1.0 : (sl - t[k]);
      exp_bg[k] = std::exp(beta * gv);
    }

    double denom = 0.0;
    for (int i = 0; i < n; i++) {
      if (Y[i] < sl)
        continue;

      double hr_ic;
      if (g_type == 0) {
        // Indicator: O(1) prefix-sum split
        double sum_q_g1 = cumq_ic[i * (m1 + 1) + cl];
        hr_ic = sum_q_g1 * exp_beta + (cumq_ic[i * (m1 + 1) + m1] - sum_q_g1);
      } else {
        // ReLU: g=0 part (k>=cl) via prefix sum; g>0 part (k<cl) loops
        hr_ic = cumq_ic[i * (m1 + 1) + m1] - cumq_ic[i * (m1 + 1) + cl];
        for (int k = 0; k < cl; k++) {
          double q_ik =
              cumq_ic[i * (m1 + 1) + k + 1] - cumq_ic[i * (m1 + 1) + k];
          hr_ic += q_ik * exp_bg[k];
        }
      }
      denom += exp_alphaX[i] * (hr_ic + q(i, m1));
    }

    if (denom > 0.0)
      h[l] = numer / denom;
  }

  return h;
}

// ============================================================
// 4. update_alphabeta_cpp
// Returns list(alpha, beta)
// ============================================================
// [[Rcpp::export]]
List update_alphabeta_cpp(NumericMatrix X,     // n x p
                          NumericVector Y,     // n
                          IntegerVector delta, // n
                          NumericVector t,     // m1 (sorted)
                          NumericVector s,     // m2 (sorted)
                          NumericMatrix q,     // n x (m1+1)
                          NumericVector alpha, // p
                          double beta,         // scalar
                          int g_type = 0       // 0 = indicator, 1 = ReLU
) {
  int n = X.nrow(), p = X.ncol();
  int m1 = t.size(), m2 = s.size();

  std::vector<double> xa(n, 0.0);
  matvec(X, alpha, xa);
  std::vector<double> exp_alphaX(n);
  for (int i = 0; i < n; i++)
    exp_alphaX[i] = std::exp(xa[i]);

  // Row-wise prefix sums of q (IC columns)
  std::vector<double> cumq_ic((size_t)n * (m1 + 1), 0.0);
  for (int i = 0; i < n; i++)
    for (int k = 0; k < m1; k++)
      cumq_ic[i * (m1 + 1) + k + 1] = cumq_ic[i * (m1 + 1) + k] + q(i, k);

  // ---- Term 1 (event contributions) ----
  std::vector<double> term1_alpha(p, 0.0);
  double term1_beta = 0.0;
  for (int i = 0; i < n; i++) {
    if (delta[i] != 1)
      continue;
    for (int j = 0; j < p; j++)
      term1_alpha[j] += X(i, j);
    // nt = number of t[k] <= Y[i]; only these have g > 0
    int nt = (int)(std::upper_bound(t.begin(), t.end(), Y[i]) - t.begin());
    double tb = 0.0;
    if (g_type == 0) {
      tb = cumq_ic[i * (m1 + 1) + nt]; // indicator: sum q[i,k] for k < nt (g=1)
    } else {
      for (int k = 0; k < nt; k++) { // relu: only k < nt contribute
        double q_ik = cumq_ic[i * (m1 + 1) + k + 1] - cumq_ic[i * (m1 + 1) + k];
        tb += q_ik * (Y[i] - t[k]);
      }
    }
    term1_beta += tb;
  }

  // ---- Term 2 and Hessian: loop over event times s[l] ----
  std::vector<double> term2_alpha(p, 0.0);
  double term2_beta = 0.0;

  arma::mat H_aa(p, p, arma::fill::zeros);
  arma::vec H_ab(p, arma::fill::zeros);
  double H_bb = 0.0;

  double exp_beta = std::exp(beta);

  // Pre-compute cl[l] = number of t[k] <= s[l] (g > 0 only for k < cl)
  std::vector<int> n_left_t(m2);
  for (int l = 0; l < m2; l++)
    n_left_t[l] = (int)(std::upper_bound(t.begin(), t.end(), s[l]) - t.begin());

  for (int l = 0; l < m2; l++) {
    double sl = s[l];

    double n_l = 0.0;
    for (int i = 0; i < n; i++)
      if (delta[i] == 1 && Y[i] == sl)
        n_l += 1.0;
    if (n_l == 0.0)
      continue;

    int cl = n_left_t[l];

    // Pre-compute g(sl, t[k]) and exp(beta*g) only for k < cl (g > 0 side)
    std::vector<double> g_vals(cl), exp_bg(cl);
    for (int k = 0; k < cl; k++) {
      double gv = (g_type == 0) ? 1.0 : (sl - t[k]);
      g_vals[k] = gv;
      exp_bg[k] = std::exp(beta * gv);
    }

    double D_l = 0.0, R1_l = 0.0, R2_l = 0.0;
    std::vector<double> S1_l(p, 0.0), XW1_l(p, 0.0);
    arma::mat S2_l(p, p, arma::fill::zeros);

    for (int i = 0; i < n; i++) {
      if (Y[i] < sl)
        continue;

      double W_il, W_il1, W_il2;
      if (g_type == 0) {
        // indicator: O(1) split — k < cl has g=1 (exp_beta), k >= cl has g=0
        // (1)
        double q_g1 = cumq_ic[i * (m1 + 1) + cl];
        double q_g0 = cumq_ic[i * (m1 + 1) + m1] - q_g1;
        W_il = q_g1 * exp_beta + q_g0;
        W_il1 = q_g1 * exp_beta; // g=1 side; g=0 side contributes 0
        W_il2 = q_g1 * exp_beta; // g^2 = g for indicator
      } else {
        // relu: k >= cl has g=0 (O(1) prefix sum); k < cl loops
        double q_g0 = cumq_ic[i * (m1 + 1) + m1] - cumq_ic[i * (m1 + 1) + cl];
        W_il = q_g0;
        W_il1 = 0.0;
        W_il2 = 0.0;
        for (int k = 0; k < cl; k++) {
          double q_ik =
              cumq_ic[i * (m1 + 1) + k + 1] - cumq_ic[i * (m1 + 1) + k];
          double gv = g_vals[k], eg = exp_bg[k];
          W_il += q_ik * eg;
          W_il1 += q_ik * gv * eg;
          W_il2 += q_ik * gv * gv * eg;
        }
      }
      W_il += q(i, m1); // RC column: g = 0, exp(beta*0) = 1

      double tw = exp_alphaX[i] * W_il;
      double tw1 = exp_alphaX[i] * W_il1;

      D_l += tw;
      R1_l += tw1;
      R2_l += exp_alphaX[i] * W_il2;
      for (int j = 0; j < p; j++) {
        S1_l[j] += X(i, j) * tw;
        XW1_l[j] += X(i, j) * tw1;
        for (int j2 = 0; j2 < p; j2++)
          S2_l(j, j2) += X(i, j) * X(i, j2) * tw;
      }
    }

    if (D_l == 0.0)
      continue;

    double rbar_l = R1_l / D_l;
    double rbar2_l = R2_l / D_l;
    arma::vec Xbar_l(p), XW1bar_l(p);
    for (int j = 0; j < p; j++) {
      Xbar_l[j] = S1_l[j] / D_l;
      XW1bar_l[j] = XW1_l[j] / D_l;
    }

    for (int j = 0; j < p; j++)
      term2_alpha[j] += n_l * Xbar_l[j];
    term2_beta += n_l * rbar_l;

    H_aa -= n_l * (S2_l / D_l - Xbar_l * Xbar_l.t());
    H_ab -= n_l * (XW1bar_l - Xbar_l * rbar_l);
    H_bb -= n_l * (rbar2_l - rbar_l * rbar_l);
  }

  // ---- Assemble score and Hessian ----
  arma::vec score(p + 1);
  for (int j = 0; j < p; j++)
    score[j] = term1_alpha[j] - term2_alpha[j];
  score[p] = term1_beta - term2_beta;

  arma::mat hessian(p + 1, p + 1, arma::fill::zeros);
  hessian.submat(0, 0, p - 1, p - 1) = H_aa;
  for (int j = 0; j < p; j++) {
    hessian(j, p) = H_ab[j];
    hessian(p, j) = H_ab[j];
  }
  hessian(p, p) = H_bb;

  arma::vec theta(p + 1);
  for (int j = 0; j < p; j++)
    theta[j] = alpha[j];
  theta[p] = beta;

  arma::vec dv;
  bool ok = arma::solve(dv, hessian, score, arma::solve_opts::no_approx);
  if (!ok)
    dv = arma::pinv(hessian) * score;
  arma::vec theta_new = theta - dv;

  NumericVector alpha_new(p);
  for (int j = 0; j < p; j++)
    alpha_new[j] = theta_new[j];

  return List::create(_["alpha"] = alpha_new, _["beta"] = theta_new[p]);
}

// ============================================================
// 5. update_gamma_cpp
// Returns p-vector gamma
// ============================================================
// [[Rcpp::export]]
NumericVector
update_gamma_cpp(NumericMatrix X,    // n x p
                 NumericVector R,    // n (right IC interval endpoint)
                 NumericVector t,    // m1 (IC time grid)
                 NumericMatrix w,    // n x m1
                 NumericVector gamma // p
) {
  int n = X.nrow(), p = X.ncol();
  int m1 = t.size();

  arma::vec score_v(p, arma::fill::zeros);
  arma::mat hessian_v(p, p, arma::fill::zeros);

  for (int k = 0; k < m1; k++) {
    double tk = t[k];

    std::vector<int> ar;
    ar.reserve(n);
    for (int i = 0; i < n; i++)
      if (R[i] >= tk)
        ar.push_back(i);
    int nr = ar.size();
    if (nr == 0)
      continue;

    std::vector<double> eGX(nr);
    for (int ii = 0; ii < nr; ii++) {
      double sg = 0.0;
      for (int j = 0; j < p; j++)
        sg += gamma[j] * X(ar[ii], j);
      eGX[ii] = std::exp(sg);
    }

    double s0 = 0.0;
    for (int ii = 0; ii < nr; ii++)
      s0 += eGX[ii];
    if (s0 == 0.0)
      continue;

    arma::vec s1(p, arma::fill::zeros);
    for (int ii = 0; ii < nr; ii++)
      for (int j = 0; j < p; j++)
        s1[j] += X(ar[ii], j) * eGX[ii];
    s1 /= s0;

    double sum_w = 0.0;
    for (int ii = 0; ii < nr; ii++)
      sum_w += w(ar[ii], k);

    for (int ii = 0; ii < nr; ii++) {
      double wik = w(ar[ii], k);
      for (int j = 0; j < p; j++)
        score_v[j] += wik * (X(ar[ii], j) - s1[j]);
    }

    arma::mat s2(p, p, arma::fill::zeros);
    for (int ii = 0; ii < nr; ii++)
      for (int j1 = 0; j1 < p; j1++)
        for (int j2 = 0; j2 < p; j2++)
          s2(j1, j2) += X(ar[ii], j1) * X(ar[ii], j2) * eGX[ii];
    s2 /= s0;

    arma::mat V = s2 - s1 * s1.t();
    hessian_v -= sum_w * V;
  }

  arma::vec gamma_v(p);
  for (int j = 0; j < p; j++)
    gamma_v[j] = gamma[j];

  arma::vec dv;
  bool ok_g = arma::solve(dv, hessian_v, score_v, arma::solve_opts::no_approx);
  if (!ok_g)
    dv = arma::pinv(hessian_v) * score_v;
  arma::vec gamma_new = gamma_v - dv;

  NumericVector result(p);
  for (int j = 0; j < p; j++)
    result[j] = gamma_new[j];
  return result;
}
