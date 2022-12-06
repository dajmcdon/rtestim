#ifndef __TF_H
#define __TF_H

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

void tf_dp(int n, double* y, double lam, double* beta);
void tf_dp_weight(int n, double* y, double* w, double lam, double* beta);

#endif
