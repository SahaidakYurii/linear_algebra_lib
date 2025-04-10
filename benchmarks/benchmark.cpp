#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "matrix.h"

using namespace std::chrono;

// Generic benchmarking function (in milliseconds)
template <typename Func>
double benchmark(Func func, int repetitions = 10) {
    double total_time = 0.0;
    for (int i = 0; i < repetitions; ++i) {
        auto start = high_resolution_clock::now();
        func();
        auto end = high_resolution_clock::now();
        total_time += duration<double, std::milli>(end - start).count();
    }
    return total_time / repetitions;
}

void run_addition_benchmark(std::ofstream& file, int size, int repetitions) {
    matrix<float> my_a(size, size), my_b(size, size);
    for (size_t i = 0; i < my_a.get_total_size(); ++i) {
        my_a.get_data()[i] = 1.0f;
        my_b.get_data()[i] = 2.0f;
    }

    Eigen::MatrixXf eigen_a = Eigen::MatrixXf::Constant(size, size, 1.0f);
    Eigen::MatrixXf eigen_b = Eigen::MatrixXf::Constant(size, size, 2.0f);

    double mylib_time = benchmark([&]() {
        matrix<float> c(my_a);
        c += my_b;
    }, repetitions);

    double eigen_time = benchmark([&]() {
        Eigen::MatrixXf c = eigen_a + eigen_b;
    }, repetitions);

    file << size << "," << mylib_time << "," << eigen_time << "\n";
}

void run_multiplication_benchmark(std::ofstream& file, int size, int repetitions) {
    matrix<float> my_a(size, size), my_b(size, size);
    for (size_t i = 0; i < my_a.get_total_size(); ++i) {
        my_a.get_data()[i] = 1.0f;
        my_b.get_data()[i] = 1.0f;
    }

    Eigen::MatrixXf eigen_a = Eigen::MatrixXf::Constant(size, size, 1.0f);
    Eigen::MatrixXf eigen_b = Eigen::MatrixXf::Constant(size, size, 1.0f);

    double mylib_time = benchmark([&]() {
        matrix<float> c = my_a * my_b;
    }, repetitions);

    double eigen_time = benchmark([&]() {
        Eigen::MatrixXf c = eigen_a * eigen_b;
    }, repetitions);

    file << size << "," << mylib_time << "," << eigen_time << "\n";
}

int main() {
    std::ofstream add_file("addition_results.csv");
    std::ofstream mult_file("multiplication_results.csv");

    add_file << "Size,YourMatrix(ms),Eigen(ms)\n";
    mult_file << "Size,YourMatrix(ms),Eigen(ms)\n";

    int repetitions = 10;

    for (int size = 100; size <= 1000; size += 100) {
        std::cout << "Running benchmarks for size: " << size << "...\n";
        run_addition_benchmark(add_file, size, repetitions);
        run_multiplication_benchmark(mult_file, size, repetitions);
    }

    add_file.close();
    mult_file.close();

    std::cout << "Benchmarks completed. Results saved to CSV files.\n";
    return 0;
}
