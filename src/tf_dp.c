/****************************************************************************
 * Copyright (C) 2014 by Taylor Arnold, Veeranjaneyulu Sadhanala,           *
 *                       Ryan Tibshirani                                    *
 *                                                                          *
 * This file is part of the glmgen library / package.                       *
 *                                                                          *
 *   glmgen is free software: you can redistribute it and/or modify it      *
 *   under the terms of the GNU Lesser General Public License as published  *
 *   by the Free Software Foundation, either version 2 of the License, or   *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   glmgen is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU Lesser General Public License for more details.                    *
 *                                                                          *
 *   You should have received a copy of the GNU Lesser General Public       *
 *   License along with glmgen. If not, see <http://www.gnu.org/licenses/>. *
 ****************************************************************************/

/**
 * @file tf_dp.c
 * @author Taylor Arnold, Veeranjaneyulu Sadhanala, Ryan Tibshirani
 * @date 2014-12-23
 * @brief Dynamic programming algorithm for the 1d fused lasso.
 *
 * Here.
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "tf_dp.h"

/**
 * @brief Implementation of Nick Johnson's dynamic programming algorithm
 * for exact O(n) calculation of the 1d fused lasso solution (at a given
 * tuning parameter value).
 * @param n                    number of observations
 * @param y                    response vector
 * @param lam                  the maximum lambda of the path
 * @param beta                 allocated space for the output
 * @return  void
 * @see tf_dp_weight
 */
void tf_dp(unsigned int n, double* y, double lam, double* beta) {
  int i;
  int k;
  int l;
  int r;
  int lo;
  int hi;
  double afirst;
  double alast;
  double bfirst;
  double blast;
  double alo;
  double blo;
  double ahi;
  double bhi;
  double* x;
  double* a;
  double* b;
  double* tm;
  double* tp;

  /* Take care of a few trivial cases */
  if (n == 0)
    return;
  if (n == 1 || lam == 0) {
    for (i = 0; i < n; i++)
      beta[i] = y[i];
    return;
  }

  x = (double*)malloc(2 * n * sizeof(double));
  a = (double*)malloc(2 * n * sizeof(double));
  b = (double*)malloc(2 * n * sizeof(double));

  /* These are the knots of the back-pointers */
  tm = (double*)malloc((n - 1) * sizeof(double));
  tp = (double*)malloc((n - 1) * sizeof(double));

  /* We step through the first iteration manually */
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
  bfirst = -y[1] - lam;
  alast = -1;
  blast = y[1] - lam;

  /* Now iterations 2 through n-1 */
  for (k = 1; k < n - 1; k++) {
    /* Compute lo: step up from l until the
       derivative is greater than -lam */
    alo = afirst;
    blo = bfirst;
    for (lo = l; lo <= r; lo++) {
      if (alo * x[lo] + blo > -lam)
        break;
      alo += a[lo];
      blo += b[lo];
    }

    /* Compute hi: step down from r until the
       derivative is less than lam */
    ahi = alast;
    bhi = blast;
    for (hi = r; hi >= lo; hi--) {
      if (-ahi * x[hi] - bhi < lam)
        break;
      ahi += a[hi];
      bhi += b[hi];
    }

    /* Compute the negative knot */
    tm[k] = (-lam - blo) / alo;
    l = lo - 1;
    x[l] = tm[k];

    /* Compute the positive knot */
    tp[k] = (lam + bhi) / (-ahi);
    r = hi + 1;
    x[r] = tp[k];

    /* Update a and b */
    a[l] = alo;
    b[l] = blo + lam;
    a[r] = ahi;
    b[r] = bhi + lam;
    afirst = 1;
    bfirst = -y[k + 1] - lam;
    alast = -1;
    blast = y[k + 1] - lam;
  }

  /* Compute the last coefficient: this is where
     the function has zero derivative */
  alo = afirst;
  blo = bfirst;
  for (lo = l; lo <= r; lo++) {
    if (alo * x[lo] + blo > 0)
      break;
    alo += a[lo];
    blo += b[lo];
  }
  beta[n - 1] = -blo / alo;

  /* Compute the rest of the coefficients, by the
     back-pointers */
  for (k = n - 2; k >= 0; k--) {
    if (beta[k + 1] > tp[k])
      beta[k] = tp[k];
    else if (beta[k + 1] < tm[k])
      beta[k] = tm[k];
    else
      beta[k] = beta[k + 1];
  }

  /* Done! Free up memory */
  free(x);
  free(a);
  free(b);
  free(tm);
  free(tp);
}

/**
 * @brief Weighted variant of the dynamic programming algorithm for the 1d
 * fused lasso problem.
 * @param n                    number of observations
 * @param y                    response vector
 * @param w                    vector of signal locations
 * @param lam                  the maximum lambda of the path
 * @param beta                 allocated space for the output
 * @return  void
 * @see tf_dp
 */
void tf_dp_weight(unsigned int n, double* y, double* w, double lam, double* beta) {
  int i;
  int k;
  int l;
  int r;
  int lo;
  int hi;
  double afirst;
  double alast;
  double bfirst;
  double blast;
  double alo;
  double blo;
  double ahi;
  double bhi;
  double* x;
  double* a;
  double* b;
  double* tm;
  double* tp;

  /* Take care of a few trivial cases */
  if (n == 0)
    return;
  if (n == 1 || lam == 0) {
    for (i = 0; i < n; i++)
      beta[i] = y[i];
    return;
  }

  /* Now deal with zero weights  */

  x = (double*)malloc(2 * n * sizeof(double));
  a = (double*)malloc(2 * n * sizeof(double));
  b = (double*)malloc(2 * n * sizeof(double));

  /* These are the knots of the back-pointers */
  tm = (double*)malloc((n - 1) * sizeof(double));
  tp = (double*)malloc((n - 1) * sizeof(double));

  /* We step through the first iteration manually */
  tm[0] = -lam / w[0] + y[0];
  tp[0] = lam / w[0] + y[0];
  l = n - 1;
  r = n;
  x[l] = tm[0];
  x[r] = tp[0];
  a[l] = w[0];
  b[l] = -w[0] * y[0] + lam;
  a[r] = -w[0];
  b[r] = w[0] * y[0] + lam;
  afirst = w[1];
  bfirst = -w[1] * y[1] - lam;
  alast = -w[1];
  blast = w[1] * y[1] - lam;

  /* Now iterations 2 through n-1 */
  for (k = 1; k < n - 1; k++) {
    /* Compute lo: step up from l until the
       derivative is greater than -lam */
    alo = afirst;
    blo = bfirst;
    for (lo = l; lo <= r; lo++) {
      if (alo * x[lo] + blo > -lam)
        break;
      alo += a[lo];
      blo += b[lo];
    }

    /* Compute hi: step down from r until the
       derivative is less than lam */
    ahi = alast;
    bhi = blast;
    for (hi = r; hi >= lo; hi--) {
      if (-ahi * x[hi] - bhi < lam)
        break;
      ahi += a[hi];
      bhi += b[hi];
    }

    /* Compute the negative knot */
    tm[k] = (-lam - blo) / alo;
    l = lo - 1;
    x[l] = tm[k];

    /* Compute the positive knot */
    tp[k] = (lam + bhi) / (-ahi);
    r = hi + 1;
    x[r] = tp[k];

    /* Update a and b */
    a[l] = alo;
    b[l] = blo + lam;
    a[r] = ahi;
    b[r] = bhi + lam;
    afirst = w[k + 1];
    bfirst = -w[k + 1] * y[k + 1] - lam;
    alast = -w[k + 1];
    blast = w[k + 1] * y[k + 1] - lam;
  }

  /* Compute the last coefficient: this is where
     the function has zero derivative */
  alo = afirst;
  blo = bfirst;
  for (lo = l; lo <= r; lo++) {
    if (alo * x[lo] + blo > 0)
      break;
    alo += a[lo];
    blo += b[lo];
  }
  beta[n - 1] = -blo / alo;

  /* Compute the rest of the coefficients, by the
     back-pointers */
  for (k = n - 2; k >= 0; k--) {
    if (beta[k + 1] > tp[k])
      beta[k] = tp[k];
    else if (beta[k + 1] < tm[k])
      beta[k] = tm[k];
    else
      beta[k] = beta[k + 1];
  }

  /* Done! Free up memory */
  free(x);
  free(a);
  free(b);
  free(tm);
  free(tp);
}

/**
 * @brief Implementation of Nick Johnson's dynamic programming algorithm
 * for exact O(n) calculation of the 1d fused lasso solution (at a given
 * tuning parameter value).
 * @param n                    number of observations
 * @param y                    response vector
 * @param past                 vector of weighted past
 * @param lam                  the maximum lambda of the path
 * @param beta                 allocated space for the output
 * @return  void
 * @see tf_dp
 */
void tf_dp_past(unsigned int n, double* y, double* past, double lam, double* beta) {
  int i;
  int k;
  int l;
  int r;
  int lo;
  int hi;
  double afirst;
  double alast;
  double bfirst;
  double blast;
  double alo;
  double blo;
  double ahi;
  double bhi;
  double* x;
  double* a;
  double* b;
  double* tm;
  double* tp;
  double* lamv;
  double* yw;
  yw = (double*)malloc(n * sizeof(double));

  for (i = 0; i < n; i++) {
    yw[i] = y[i] / past[i];
  }
  /* Take care of a few trivial cases */
  if (n == 0)
    return;
  if (n == 1 || lam == 0) {
    for (i = 0; i < n; i++)
      beta[i] = yw[i];
    return;
  }

  x = (double*)malloc(2 * n * sizeof(double));
  a = (double*)malloc(2 * n * sizeof(double));
  b = (double*)malloc(2 * n * sizeof(double));

  /* These are the knots of the back-pointers */
  tm = (double*)malloc((n - 1) * sizeof(double));
  tp = (double*)malloc((n - 1) * sizeof(double));

  lamv = (double*)malloc(n * sizeof(double));
  for (i = 0; i < n; i++) {
    lamv[i] = lam / past[i];
  }

  /* We step through the first iteration manually */
  tm[0] = -lamv[0] + yw[0];
  tp[0] = lamv[0] + yw[0];
  l = n - 1;
  r = n;
  x[l] = tm[0];
  x[r] = tp[0];
  a[l] = 1;
  b[l] = -yw[0] + lamv[1];
  a[r] = -1;
  b[r] = yw[0] + lamv[1];
  afirst = 1;
  bfirst = -yw[1] - lamv[1];
  alast = -1;
  blast = yw[1] - lamv[1];

  /* Now iterations 2 through n-1 */
  for (k = 1; k < n - 1; k++) {
    /* Compute lo: step up from l until the
       derivative is greater than -lam */
    alo = afirst;
    blo = bfirst;
    for (lo = l; lo <= r; lo++) {
      if (alo * x[lo] + blo > -lamv[k + 1])
        break;
      alo += a[lo];
      blo += b[lo];
    }

    /* Compute hi: step down from r until the
       derivative is less than lam */
    ahi = alast;
    bhi = blast;
    for (hi = r; hi >= lo; hi--) {
      if (-ahi * x[hi] - bhi < lamv[k + 1])
        break;
      ahi += a[hi];
      bhi += b[hi];
    }

    /* Compute the negative knot */
    tm[k] = (-lamv[k + 1] - blo) / alo;
    l = lo - 1;
    x[l] = tm[k];

    /* Compute the positive knot */
    tp[k] = (lamv[k + 1] + bhi) / (-ahi);
    r = hi + 1;
    x[r] = tp[k];

    /* Update a and b */
    a[l] = alo;
    b[l] = blo + lamv[k + 1];
    a[r] = ahi;
    b[r] = bhi + lamv[k + 1];
    afirst = 1;
    bfirst = -yw[k + 1] - lamv[k + 1];
    alast = -1;
    blast = yw[k + 1] - lamv[k + 1];
  }

  /* Compute the last coefficient: this is where
     the function has zero derivative */
  alo = afirst;
  blo = bfirst;
  for (lo = l; lo <= r; lo++) {
    if (alo * x[lo] + blo > 0)
      break;
    alo += a[lo];
    blo += b[lo];
  }
  beta[n - 1] = -blo / alo;

  /* Compute the rest of the coefficients, by the
     back-pointers */
  for (k = n - 2; k >= 0; k--) {
    if (beta[k + 1] > tp[k])
      beta[k] = tp[k];
    else if (beta[k + 1] < tm[k])
      beta[k] = tm[k];
    else
      beta[k] = beta[k + 1];
  }

  /* Done! Free up memory */
  free(x);
  free(a);
  free(b);
  free(tm);
  free(tp);
  free(yw);
  free(lamv);
}

/**
 * @brief Weighted variant of the dynamic programming algorithm for the 1d
 * fused lasso problem. This function is modified to handle weighted signals.
 * @param n                    number of observations
 * @param y                    response vector
 * @param past                 vector of weighted past
 * @param w                    vector of signal locations
 * @param lam                  the maximum lambda of the path
 * @param beta                 allocated space for the output
 * @return  void
 * @see tf_dp_weight
 */
void tf_dp_past_weight(unsigned int n,
                       double* y,
                       double* past,
                       double* w,
                       double lam,
                       double* beta) {
  int i;
  int k;
  int l;
  int r;
  int lo;
  int hi;
  double afirst;
  double alast;
  double bfirst;
  double blast;
  double alo;
  double blo;
  double ahi;
  double bhi;
  double* x;
  double* a;
  double* b;
  double* tm;
  double* tp;
  double* lamv;

  for (i = 0; i < n; i++) {
    y[i] /= past[i];
  }
  /* Take care of a few trivial cases */
  if (n == 0)
    return;
  if (n == 1 || lam == 0) {
    for (i = 0; i < n; i++)
      beta[i] = y[i];
    return;
  }

  /* Now deal with zero weights  */

  x = (double*)malloc(2 * n * sizeof(double));
  a = (double*)malloc(2 * n * sizeof(double));
  b = (double*)malloc(2 * n * sizeof(double));

  /* These are the knots of the back-pointers */
  tm = (double*)malloc((n - 1) * sizeof(double));
  tp = (double*)malloc((n - 1) * sizeof(double));

  lamv = (double*)malloc((n) * sizeof(double));
  for (i = 0; i < n; i++) {
    lamv[i] = lam / past[i];
  }

  /* We step through the first iteration manually */
  tm[0] = -lamv[0] / w[0] + y[0];
  tp[0] = lamv[0] / w[0] + y[0];
  l = n - 1;
  r = n;
  x[l] = tm[0];
  x[r] = tp[0];
  a[l] = w[0];
  b[l] = -w[0] * y[0] + lamv[1];
  a[r] = -w[0];
  b[r] = w[0] * y[0] + lamv[1];
  afirst = w[1];
  bfirst = -w[1] * y[1] - lamv[1];
  alast = -w[1];
  blast = w[1] * y[1] - lamv[1];

  /* Now iterations 2 through n-1 */
  for (k = 1; k < n - 1; k++) {
    /* Compute lo: step up from l until the
       derivative is greater than -lam */
    alo = afirst;
    blo = bfirst;
    for (lo = l; lo <= r; lo++) {
      if (alo * x[lo] + blo > -lamv[k + 1])
        break;
      alo += a[lo];
      blo += b[lo];
    }

    /* Compute hi: step down from r until the
       derivative is less than lam */
    ahi = alast;
    bhi = blast;
    for (hi = r; hi >= lo; hi--) {
      if (-ahi * x[hi] - bhi < lamv[k + 1])
        break;
      ahi += a[hi];
      bhi += b[hi];
    }

    /* Compute the negative knot */
    tm[k] = (-lamv[k + 1] - blo) / alo;
    l = lo - 1;
    x[l] = tm[k];

    /* Compute the positive knot */
    tp[k] = (lamv[k + 1] + bhi) / (-ahi);
    r = hi + 1;
    x[r] = tp[k];

    /* Update a and b */
    a[l] = alo;
    b[l] = blo + lamv[k + 1];
    a[r] = ahi;
    b[r] = bhi + lamv[k + 1];
    afirst = w[k + 1];
    bfirst = -w[k + 1] * y[k + 1] - lamv[k + 1];
    alast = -w[k + 1];
    blast = w[k + 1] * y[k + 1] - lamv[k + 1];
  }

  /* Compute the last coefficient: this is where
     the function has zero derivative */
  alo = afirst;
  blo = bfirst;
  for (lo = l; lo <= r; lo++) {
    if (alo * x[lo] + blo > 0)
      break;
    alo += a[lo];
    blo += b[lo];
  }
  beta[n - 1] = -blo / alo;

  /* Compute the rest of the coefficients, by the
     back-pointers */
  for (k = n - 2; k >= 0; k--) {
    if (beta[k + 1] > tp[k])
      beta[k] = tp[k];
    else if (beta[k + 1] < tm[k])
      beta[k] = tm[k];
    else
      beta[k] = beta[k + 1];
  }

  /* Done! Free up memory */
  free(x);
  free(a);
  free(b);
  free(tm);
  free(tp);
}
