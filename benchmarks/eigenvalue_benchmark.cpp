#include <iostream>
#include <chrono>
#include <fstream>
#include "linalg.h"
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include "squareMatrix.h"
#include "benchmark_utils.h"

using namespace std;

linalg::squareMatrix<double> generate_symmetric_matrix(size_t n) {
    linalg::vector<double> data(n * n, 0);
    std::mt19937 gen(random_device{}());
    std::uniform_real_distribution<> dis(-10.0, 10.0);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i; j < n; ++j) {
            double val = dis(gen);
            data[i * n + j] = val;
            data[j * n + i] = val;
        }
    }

    return linalg::squareMatrix<double>(n, data);
}

void run_benchmark(const std::string& output_csv, int min_size = 5, int max_size = 75, int step = 5, int repetitions = 10) {
    ofstream file(output_csv);
    file << "Size,YourMatrix(ms),Eigen(ms)\n";

    for (int size = min_size; size <= max_size; size += step) {
        cout << "Benchmarking size: " << size << "..." << endl;

        double total_your_time = 0;
        double total_eigen_time = 0;

        linalg::squareMatrix<double> my_mat = generate_symmetric_matrix(size);

        Eigen::MatrixXd eigen_mat(size, size);
        for (int r = 0; r < size; ++r)
            for (int c = 0; c < size; ++c)
                eigen_mat(r, c) = my_mat(r, c);

        total_your_time += benchmark([&]() {
            volatile auto eigvals = my_mat.eigenvalues();
        }, 1);

        total_eigen_time += benchmark([&]() {
            Eigen::EigenSolver<Eigen::MatrixXd> solver(eigen_mat);
            volatile auto eigvals = solver.eigenvalues();
        }, 1);

        double avg_your_time = total_your_time / repetitions;
        double avg_eigen_time = total_eigen_time / repetitions;

        file << size << "," << avg_your_time << "," << avg_eigen_time << "\n";
    }

    file.close();
    cout << "Benchmark complete! Results saved to: " << output_csv << endl;
}

int main() {
    run_benchmark("eigenvalue_benchmark.csv");
    return 0;
}
