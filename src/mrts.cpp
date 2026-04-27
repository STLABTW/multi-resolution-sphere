// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#define _USE_MATH_DEFINES
#include <RcppEigen.h>
#include <Rcpp.h>
#include <RcppNumerical.h>
#include <cmath>
#include <vector>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Numer;
using namespace Rcpp;

// Integrand: log(1 - x) / x, with the removable singularity at x = 0.
class Func_logterm : public Func {
public:
    double operator()(const double& x) const {
        if (x == 0.0) return -1.0;
        return std::log(1.0 - x) / x;
    }
};

// Kernel value K(p1, p2) for two points on the sphere given in
// (latitude, longitude) degrees. Used to build the n x n knot kernel.
// [[Rcpp::export]]
double cpp_Kf(double L1, double l1, double L2, double l2) {
    const double mia = M_PI / 180.0;
    double a = std::sin(L1 * mia) * std::sin(L2 * mia) +
               std::cos(L1 * mia) * std::cos(L2 * mia) *
                   std::cos(l1 * mia - l2 * mia);
    double b = std::max(-1.0, std::min(a, 1.0));
    double aaa = std::acos(b);

    if (std::cos(aaa) == -1.0) {
        return 1.0 - std::pow(M_PI, 2) / 6.0;
    }
    double aa = 0.5 + std::cos(aaa) / 2.0;
    double lower = 0.0;
    Func_logterm f;
    double err_est;
    int err_code;
    double res = integrate(f, lower, aa, err_est, err_code);
    return 1.0 - std::pow(M_PI, 2) / 6.0 - res;
}

// n x n kernel matrix on the knot set given as latitude/longitude vectors.
// [[Rcpp::export]]
NumericMatrix cpp_K(NumericVector X, NumericVector Y, int n) {
    NumericMatrix xx(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            xx(i, j) = cpp_Kf(X[i], Y[i], X[j], Y[j]);
        }
    }
    return xx;
}

// Linear interpolation lookup against a sorted grid.
static double interp1(double aa,
                      const NumericVector& aa_grid,
                      const NumericVector& tab_vals) {
    int n = aa_grid.size();
    if (aa <= aa_grid[0]) return tab_vals[0];
    if (aa >= aa_grid[n - 1]) return tab_vals[n - 1];

    int left = 0, right = n - 1;
    while (right - left > 1) {
        int mid = (left + right) / 2;
        if (aa_grid[mid] <= aa) left = mid;
        else right = mid;
    }
    double t = (aa - aa_grid[left]) / (aa_grid[left + 1] - aa_grid[left]);
    return tab_vals[left] * (1.0 - t) + tab_vals[left + 1] * t;
}

// MRTS basis matrix evaluated on the prediction grid.
// Uses precomputed integral lookup table (aa_grid, tab_vals) and the
// eigen-projection matrix eiKvecmval = eig.vec / eig.val.
// [[Rcpp::export]]
NumericMatrix cpp_Kmatrix6(int KK,
                           NumericMatrix X,
                           NumericMatrix ggrids,
                           NumericVector Konev,
                           NumericMatrix eiKvecmval,
                           int n,
                           int N,
                           NumericVector aa_grid,
                           NumericVector tab_vals) {
    NumericMatrix xx(N, KK);
    const double mia = M_PI / 180.0;

    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (int i = 0; i < N; i++) {
        double L1 = ggrids(i, 0);
        double l1 = ggrids(i, 1);

        std::vector<double> f2(n);
        for (int j = 0; j < n; j++) {
            double L2 = X(j, 0);
            double l2 = X(j, 1);

            double a = std::sin(L1 * mia) * std::sin(L2 * mia) +
                       std::cos(L1 * mia) * std::cos(L2 * mia) *
                           std::cos(l1 * mia - l2 * mia);
            double b = std::max(-1.0, std::min(a, 1.0));
            double aaa = std::acos(b);
            double result;

            if (std::cos(aaa) == -1.0) {
                result = 1.0 - std::pow(M_PI, 2) / 6.0;
            } else {
                double aa = 0.5 + std::cos(aaa) / 2.0;
                double res = interp1(aa, aa_grid, tab_vals);
                result = 1.0 - std::pow(M_PI, 2) / 6.0 - res;
            }
            f2[j] = result;
        }

        std::vector<double> t(n);
        for (int j = 0; j < n; j++) {
            t[j] = f2[j] - Konev[j];
        }

        xx(i, 0) = std::sqrt(1.0 / static_cast<double>(n));
        for (int k = 1; k < KK; k++) {
            double s = 0.0;
            for (int j = 0; j < n; j++) {
                s += t[j] * eiKvecmval(j, k - 1);
            }
            xx(i, k) = s;
        }
    }
    return xx;
}
