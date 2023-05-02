#ifndef __TF_DP_H
#define __TF_DP_H

/* Dynamic programming routines */
void tf_dp (int n, double *y, double lam, double *beta);
void tf_dp_weight (int n, double *y, double *w, double lam, double *beta);
void tf_dp_past(int n, double* y, double* past, double lam, double* beta);
void tf_dp_past_weight(int n, double* y, double* past, double* w, double lam,
                       double* beta);

#endif
