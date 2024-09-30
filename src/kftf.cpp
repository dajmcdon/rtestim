#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <dspline.h>
#include "utils.h"
#include "kftf.h"

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(dspline)]]

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
using Eigen::RowVectorXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Rcpp::_;
using Rcpp::IntegerVector;
using Rcpp::List;
using Rcpp::Named;
using Rcpp::NumericVector;

// Always prefix the types in function signatures
void f1step(double y,
            double c,
            double Z,
            double H,
            const Eigen::MatrixXd& A,
            double RQR,
            Eigen::VectorXd& a,
            Eigen::MatrixXd& P,
            double& vt,
            double& Ft,
            Eigen::VectorXd& Kt) {
  VectorXd a_temp = A * a;
  a_temp(0) += c;
  MatrixXd Ptemp = computePtemp(A, P);  // A * P * A.transpose();
  Ptemp(0, 0) += RQR;

  vt = y - Z * a_temp(0);
  Ft = pow(Z, 2) * Ptemp(0, 0) + H;
  Kt = Ptemp.col(0) * Z;
  a = a_temp + Kt * vt / Ft;
  P = Ptemp - Kt * Z * Ptemp.row(0) / Ft;

  // symmetrize
  Ptemp = P;
  Ptemp += P.transpose();
  P = Ptemp / 2;

  // Some entries of P can be _really_ small in magnitude, set them to 0
  // but maintain symmetric / posdef
  for (int i = 0; i < P.rows(); i++) {
    for (int j = 0; j <= i; j++) {
      if (abs(P(i, j)) < 1e-30) {
        if (i == j) {
          P.row(i).setZero();
          P.col(i).setZero();
        } else {
          P(i, j) = 0;
          P(j, i) = 0;
        }
      }
    }
  }
}

void df1step(double y,
             double Z,
             double H,
             const Eigen::MatrixXd& A,
             double RQR,
             Eigen::VectorXd& a,
             Eigen::MatrixXd& P,
             Eigen::MatrixXd& Pinf,
             int& rankp,
             double& vt,
             double& Ft,
             double& Finf,
             Eigen::VectorXd& Kt,
             Eigen::VectorXd& Kinf) {
  double tol = Eigen::NumTraits<double>::epsilon();
  tol = std::sqrt(tol);
  int k = a.size();
  MatrixXd Ptemp(k, k);

  a = A * a;
  P = computePtemp(A, P);  // A * P * A.transpose();
  P(0, 0) += RQR;
  Pinf = computePtemp(A, Pinf);  // A * Pinf * A.transpose();

  vt = y - Z * a(0);
  Kt = P.col(0) * Z;
  Ft = Z * Kt(0) + H;
  Kinf = Pinf.col(0) * Z;
  Finf = Z * Kinf(0);

  if (Finf > tol) {  // should always happen
    a += vt * Kinf / Finf;
    P += Ft * Kinf * Kinf.transpose() / pow(Finf, 2);
    P -= Kt * Kinf.transpose() / Finf + Kinf * Kt.transpose() / Finf;
    Pinf -= Kinf * Kinf.transpose() / Finf;
    rankp--;
  } else {  // should never happen
    Finf = 0;
    if (Ft > tol) {
      a += vt * Kt / Ft;
      P -= Kt * Kt.transpose() / Ft;
    }
  }
  if (Ft < tol)
    Ft = 0;

  // symmetrize
  Ptemp = P;
  Ptemp += P.transpose();
  P = Ptemp / 2;
  Ptemp = Pinf;
  Ptemp += Pinf.transpose();
  Pinf = Ptemp / 2;

  // Fix possible negative definiteness, should never happen
  for (int i = 0; i < k; i++) {
    if (P(i, i) < 0) {
      P.row(i).setZero();
      P.col(i).setZero();
    }
    if (Pinf(i, i) < 0) {
      Pinf.row(i).setZero();
      Pinf.col(i).setZero();
    }
  }
}

void kftf(const Eigen::VectorXd& y,
          int k,
          const Eigen::VectorXd& weights,
          const Eigen::VectorXd& x,
          const Eigen::VectorXd& c,
          double lambda,
          Eigen::VectorXd& theta) {
  VectorXd wsqrt = weights.array().sqrt();
  VectorXd yw = y.array() * wsqrt.array();
  int n = y.size();

  // define transition matrix A
  MatrixXd A = MatrixXd::Zero(k, k);
  A.block(1, 0, k - 1, k - 1).diagonal() = VectorXd::Ones(k - 1);
  //double s = -1;

  // construct dense matrix Dseq 
  // MatrixXd Dseq = b_mat(k, x, VectorXi::LinSpaced(n - k, 0, n - k - 1));
  // or, convert sparse matrix to dense matrix Dseq 
  NumericVector xn(Rcpp::wrap(x));
  Eigen::SparseMatrix Dseq_smat = dspline::rcpp_b_mat(k, xn, false, Rcpp::seq(k, n - 1), false);
  MatrixXd Dseq = smat_to_mat(Dseq_smat, k,  false);
  
  VectorXd RQR = VectorXd::Zero(n - k);
  VectorXd a1 = VectorXd::Zero(k);
  MatrixXd at = MatrixXd::Zero(k, n + 1);
  MatrixXd P1 = MatrixXd::Zero(k, k);
  MatrixXd Pt = MatrixXd::Zero(k * k, n + 1);
  MatrixXd P1inf = MatrixXd::Identity(k, k);
  MatrixXd Pinf = MatrixXd::Zero(k * k, n + 1);
  Pinf.col(0) = Map<Eigen::VectorXd>(P1inf.data(), k * k);

  // forward
  int d = 0;
  int rankp = k;
  double vt_b = 0.0;
  double Ft_b = 0.0;
  double Finf_b = 0.0;
  VectorXd vt = VectorXd::Zero(n);
  VectorXd Ft = VectorXd::Zero(n);
  VectorXd Finf = VectorXd::Zero(n);
  MatrixXd Kt = MatrixXd::Zero(k, n);
  MatrixXd Kinf = MatrixXd::Zero(k, n);
  VectorXd Kt_b = VectorXd::Zero(k);
  VectorXd Kinf_b = VectorXd::Zero(k);

  // re-construct Dseq in which each row saves values for A.row(0) per iterate
  VectorXd s_seq = VectorXd::Ones(n - k);
  s_seq = Dseq.block(0, k, n - k, 1);
  Dseq.conservativeResize(n - k, k);
  Eigen::MatrixXd firstRow(1, k);
  for (int i = 0; i < n - k; i++) {
    firstRow = -Dseq.row(i) / s_seq(i);
    std::reverse(firstRow.data(), firstRow.data() + k);
    Dseq.row(i) = firstRow;
  }

  RQR = s_seq.array().square().inverse();
  RQR *= 1 / lambda;

  A.row(0) = Dseq.row(0);
  while (rankp > 0 && d < n) {
    df1step(yw(d), wsqrt(d), 1, A, RQR(0), a1, P1, P1inf, rankp, vt_b, Ft_b,
            Finf_b, Kt_b, Kinf_b);
    at.col(d + 1) = a1;
    vt(d) = vt_b;
    Ft(d) = Ft_b;
    Finf(d) = Finf_b;
    Kt.col(d) = Kt_b;
    Kinf.col(d) = Kinf_b;
    Pt.col(d + 1) = Map<Eigen::VectorXd>(P1.data(), k * k);
    Pinf.col(d + 1) = Map<Eigen::VectorXd>(P1inf.data(), k * k);
    d++;
  }

  for (int i = d; i < n; i++) {
    f1step(yw(i), c(i - d) / s_seq(i - d), wsqrt(i), 1, A, RQR(i - d), a1, P1,
           vt_b, Ft_b, Kt_b);
    vt(i) = vt_b;
    Ft(i) = Ft_b;
    Kt.col(i) = Kt_b;
    at.col(i + 1) = a1;
    Pt.col(i + 1) = Map<VectorXd>(P1.data(), k * k);
    // update A for next iterate:
    if (i < n - 1)
      A.row(0) = Dseq.row(i - d + 1);
  }

  // backward
  VectorXd r = VectorXd::Zero(k);
  VectorXd r1 = VectorXd::Zero(k);
  VectorXd rtmp = VectorXd::Zero(k);
  MatrixXd L0 = MatrixXd::Zero(k, k);
  MatrixXd L1 = MatrixXd::Zero(k, k);
  MatrixXd Ptemp = MatrixXd::Zero(k, k);

  for (int i = n - 1; i >= d; i--) {
    P1 = Map<MatrixXd>(Pt.col(i + 1).data(), k, k);
    theta(i) = at.col(i + 1)(0) - P1.row(0) * r;

    L0.col(0) = Kt.col(i) * wsqrt(i) / Ft(i);
    r1 = r - L0.transpose() * r;
    r1(0) -= wsqrt(i) * vt(i) / Ft(i);
    r = A.transpose() * r1;
    // update A for next iterate:
    if (i > d)
      A.row(0) = Dseq.row(i - d - 1);
  }

  for (int i = d - 1; i >= 0; i--) {
    P1 = Map<MatrixXd>(Pt.col(i + 1).data(), k, k);
    P1inf = Map<MatrixXd>(Pinf.col(i + 1).data(), k, k);
    a1 = at.col(i + 1);
    theta(i) = a1(0) - P1.row(0) * r - P1inf.row(0) * r1;

    if (i > 0) {
      if (Finf(i) > 0) {
        // simple version w/o mat multiplication:
        L0.col(0) = Kinf.col(i) * wsqrt(i) / Finf(i);
        L1.col(0) = L0.col(0) * Ft(i) / Finf(i);
        L1.col(0) -= Kt.col(i) * wsqrt(i) / Finf(i);
        rtmp = r1 - L0.transpose() * r1 + L1.transpose() * r;
        rtmp(0) -= wsqrt(i) * vt(i) / Finf(i);
        r1 = A.transpose() * rtmp;
        rtmp = r - L0.transpose() * r;
        r = A.transpose() * rtmp;
      } else {
        L1.col(0) = Kt.col(i) * wsqrt(i) / Ft(i);
        rtmp = r - L1.transpose() * r;
        rtmp(0) -= wsqrt(i) * vt(i) / Ft(i);
        r = A.transpose() * rtmp;
        rtmp = r1 - L1.transpose() * r1;
        r1 = A.transpose() * rtmp;
      }
    }
  }
}

// [[Rcpp::export]]
Eigen::VectorXd kftf_test(Eigen::VectorXd y,
                          int k,
                          Eigen::VectorXd weights,
                          Eigen::VectorXd x,
                          Eigen::VectorXd c,
                          double lambda) {
  VectorXd theta(y.size());
  kftf(y, k, weights, x, c, lambda, theta);

  return theta;
}
