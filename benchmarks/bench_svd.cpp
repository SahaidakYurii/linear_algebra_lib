// bench_svd_only_my_lib.cpp
#include "squareMatrix.h"
#include <chrono>
#include <tuple>
#include <cmath>
#include <iostream>

using Clock    = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;
using namespace linalg;  // brings in your vector<T> and squareMatrix<T>

int main() {
    // 1) Build a simple 256×256 test matrix (identity)
    size_t n = 256;
    vector<double> data(n*n, 0.0);
    for (size_t i = 0; i < n; ++i)
        data[i*n + i] = 1.0;

    // direct‐init calls your (size, vector<T>) ctor
    squareMatrix<double> A(n, data);

    // 2) Time your SVD
    auto t0 = Clock::now();
    auto svd_res = A.svd();      // returns tuple<U,S,VT>
    auto t1 = Clock::now();
    double my_time = Duration(t1 - t0).count();

    // 3) Extract U, S, Vᵀ via your newly unambiguous copy‐ctor
    squareMatrix<double> U  = std::get<0>(svd_res);
    squareMatrix<double> S  = std::get<1>(svd_res);
    squareMatrix<double> VT = std::get<2>(svd_res);

    // 4) Reconstruct A and measure Frobenius‐error
    auto Arec = U.multiply(S).multiply(VT);
    double err = 0;
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) {
            double d = A(i,j) - Arec(i,j);
            err += d*d;
        }
    err = std::sqrt(err);

    // 5) Print results
    std::cout
            << "SVD time           = " << my_time << " s\n"
            << "Reconstruction err = " << err      << "\n";

    return 0;
}
