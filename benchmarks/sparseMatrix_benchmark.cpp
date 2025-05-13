#include <iostream>
#include <fstream>
#include <random>
#include <vector>

#include "linalg.h"
#include "benchmark_utils.h"

static linalg::vector<float> make_data(const size_t size, const float density) {
    linalg::vector<float> data(size * size, 0.0f);
    std::mt19937_64 eng(42);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    for (size_t i = 0; i < data.size(); ++i) {
        if (dist(eng) < density) data[i] = 1.0f;
    }
    return data;
}

void run_sparse_vs_dense_mult(std::ofstream& file,
                              const size_t size,
                              const float density,
                              const size_t repetitions)
{
    const auto dataA = make_data(size, density);
    const auto dataB = make_data(size, density);

    const linalg::matrix denseA(size, size, dataA);
    const linalg::matrix denseB(size, size, dataB);

    const linalg::sparseMatrix sparseA(size, size, dataA);
    const linalg::sparseMatrix sparseB(size, size, dataB);

    const double dense_t = benchmark([&]() {
        const linalg::matrix<float> C = denseA * denseB;
    }, repetitions);

    const double sparse_t = benchmark([&]() {
        const linalg::sparseMatrix<float> C = sparseA * sparseB;
    }, repetitions);

    file << size << ","
         << density << ","
         << dense_t << ","
         << sparse_t << "\n";
}

int main() {
    std::ofstream csv("sparse_vs_dense_mult.csv");
    csv << "Size,Density,DenseMult(ms),SparseMult(ms)\n";

    const int repetitions = 1;
    const std::vector<int> sizes   = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000};
    const std::vector<float> densities = {0.01f, 0.02f, 0.03f, 0.4f, 0.05f, 0.1f, 0.2f, 0.5f, 1.0f};

    for (const int sz : sizes) {
        std::cout << "Size=" << sz << "  ";
        for (const float d : densities) {
            std::cout << "ρ=" << d << " ";
            run_sparse_vs_dense_mult(csv, sz, d, repetitions);
        }
        std::cout << "\n";
    }

    csv.close();
    std::cout << "Done. Results in sparse_vs_dense_mult.csv\n";
    return 0;
}
