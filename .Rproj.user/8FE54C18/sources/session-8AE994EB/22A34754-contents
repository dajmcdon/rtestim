#include <RcppArmadillo.h>
#include <boost/math/special_functions/lambert_w.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins("cpp11")]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
vec prox_dp_norm(vec &y, double &lam)
{
    // Take care of a few trivial cases
    int n = y.size();
    if (n == 0)
        return (y);
    vec theta(n);
    if (n == 1 || lam == 0) {
        for (int i = 0; i < n; i++) theta[i] = y[i];
        return (theta);
    }

    // These are used to store the derivative of the
    // piecewise linear function of interest
    double afirst, alast, bfirst, blast;
    double *x = (double *)malloc(2 * n * sizeof(double));
    double *a = (double *)malloc(2 * n * sizeof(double));
    double *b = (double *)malloc(2 * n * sizeof(double));
    int l, r;

    // These are the knots of the back-pointers
    double *tm = (double *)malloc((n - 1) * sizeof(double));
    double *tp = (double *)malloc((n - 1) * sizeof(double));

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
            if (alo * x[lo] + blo > -lam) break;
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
            if (-ahi * x[hi] - bhi < lam) break;
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
    for (lo = l; lo <= r; lo++)
    {
        if (alo * x[lo] + blo > 0)
            break;
        alo += a[lo];
        blo += b[lo];
    }
    theta[n - 1] = -blo / alo;

    // Compute the rest of the coefficients, by the
    // back-pointers
    for (int k = n - 2; k >= 0; k--)
    {
        if (theta[k + 1] > tp[k])
            theta[k] = tp[k];
        else if (theta[k + 1] < tm[k])
            theta[k] = tm[k];
        else
            theta[k] = theta[k + 1];
    }

    return theta;
}



double update_primal(double c, double mu2) {
    c -= boost::math::lambert_w0(exp(c) / mu2);
    return c;
}

// [[Rcpp::export]]
List admm(int M, vec y, vec x, int n, vec theta, vec z,
          vec u, double lambda, double rho, double mu, sp_mat D,
          double tol = 1e-3)
{

    int iter;
    double r_norm, s_norm;
    vec z_old = z;
    sp_mat Dt = D.t();
    sp_mat DD = Dt * D;
    double mu2 = n * mu;

    // start of iteration:
    for (iter = 0; iter < M; iter++)
    {
        // update primal variable:
        vec vec_c = y / n - rho * DD * theta + rho * Dt * (z - u) + mu * theta;
        vec_c = vec_c / mu + log(x);
        vec_c.transform([&](double vec_c)
                        { return update_primal(vec_c, mu2); });
        theta = vec_c - log(x);

        // update alternating variable:
        vec y_dp = D * theta + u;
        double lam_dp = lambda / rho;
        z = prox_dp_norm(y_dp, lam_dp);

        // update dual variable:
        u += D * theta - z;

        // stopping criteria check:
        vec r = D * theta - z;
        r_norm = sqrt(mean(square(r)));
        // dual residuals:
        vec s = z_old - z;
        s_norm = sqrt(mean(square(s)));

        if ((r_norm < tol) && (s_norm < tol))
        {
            iter++;
            break;
        }
        // auxiliary variables update:
        z_old = z;
    }

    return List::create(
        Named("theta") = wrap(theta), Named("z") = wrap(z),
        Named("u") = wrap(u), Named("prim_res") = r_norm,
        Named("dual_res") = s_norm, Named("iter_num") = iter);
}
