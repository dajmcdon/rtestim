#include <RcppArmadillo.h>
#include "dp.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;

vec prox_dp_norm(vec y, double lam) {
  // Take care of a few trivial cases
  int n = y.size();
  if (n == 0)
    return (y);
  vec theta(n);
  if (n == 1 || lam == 0) {
    for (int i = 0; i < n; i++)
      theta[i] = y[i];
    return theta;
  }

  // These are used to store the derivative of the
  // piecewise linear function of interest
  double afirst, alast, bfirst, blast;
  vec x(2 * n);
  vec a(2 * n);
  vec b(2 * n);
  int l, r;

  // These are the knots of the back-pointers
  vec tm(n - 1);
  vec tp(n - 1);

  // We step through the first iteration manually
  tm[0] = -lam + y[0];
  tp[0] = lam + y[0];
  l = n - 1;
  r = n;
  x[l] = tm[0];
  x[r] = tp[0];
  a[l] = 1;
  b[l] = -y[0] + lam;
  a[r] = -1;
  b[r] = y[0] + lam;
  afirst = 1;
  bfirst = -lam - y[1];
  alast = -1;
  blast = -lam + y[1];

  // Now iterations 2 through n-1
  int lo, hi;
  double alo, blo, ahi, bhi;
  for (int k = 1; k < n - 1; k++) {
    // Compute lo: step up from l until the
    // derivative is greater than -lam
    alo = afirst;
    blo = bfirst;
    for (lo = l; lo <= r; lo++) {
      if (alo * x[lo] + blo > -lam)
        break;
      alo += a[lo];
      blo += b[lo];
    }

    // Compute the negative knot
    tm[k] = (-lam - blo) / alo;
    l = lo - 1;
    x[l] = tm[k];

    // Compute hi: step down from r until the
    // derivative is less than lam
    ahi = alast;
    bhi = blast;
    for (hi = r; hi >= l; hi--) {
      if (-ahi * x[hi] - bhi < lam)
        break;
      ahi += a[hi];
      bhi += b[hi];
    }

    // Compute the positive knot
    tp[k] = (lam + bhi) / (-ahi);
    r = hi + 1;
    x[r] = tp[k];

    // Update a and b
    a[l] = alo;
    b[l] = blo + lam;
    a[r] = ahi;
    b[r] = bhi + lam;
    afirst = 1;
    bfirst = -lam - y[k + 1];
    alast = -1;
    blast = -lam + y[k + 1];
  }

  // Compute the last coefficient: this is where
  // the function has zero derivative

  alo = afirst;
  blo = bfirst;
  for (lo = l; lo <= r; lo++) {
    if (alo * x[lo] + blo > 0)
      break;
    alo += a[lo];
    blo += b[lo];
  }
  theta[n - 1] = -blo / alo;

  // Compute the rest of the coefficients, by the
  // back-pointers
  for (int k = n - 2; k >= 0; k--) {
    if (theta[k + 1] > tp[k])
      theta[k] = tp[k];
    else if (theta[k + 1] < tm[k])
      theta[k] = tm[k];
    else
      theta[k] = theta[k + 1];
  }

  // Done!
  return theta;
}