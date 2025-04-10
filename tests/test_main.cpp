// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#include "linalg.h"

#include <cassert>
#include <cmath>
#include <iostream>


 template<typename T, size_t N>
 void print_matrix(const tensor<T, N>& matrix) {
     auto shape = matrix.shape();
     size_t index = 0;
     for (size_t i = 0; i < shape[0]; ++i) {
         for (size_t j = 0; j < shape[1]; ++j) {
             std::cout << matrix.data_m[index++] << " ";
         }
         std::cout << std::endl;
     }
 }

 template<typename T, size_t N>
 void assert_matrix_equal(const tensor<T, N>& result, const tensor<T, N>& expected, const std::string& test_name) {
     if (result.shape() != expected.shape()) {
         std::cout << "FAIL [" << test_name << "]: Dimension mismatch." << std::endl;
         return;
     }
     const auto& result_data = result.get_data();
     const auto& expected_data = expected.get_data();
     for (size_t i = 0; i < result.get_total_size(); ++i) {
         if (result_data[i] != expected_data[i]) {
             std::cout << "FAIL [" << test_name << "]: Values do not match at index " << i << ". Expected " << expected_data[i] << ", got " << result_data[i] << "." << std::endl;
             return;
         }
     }
     std::cout << "PASS [" << test_name << "]" << std::endl;
 }


 void test_addition() {
     tensor<int, 2> mat1({1, 2, 3, 4}, 2, 2);
     tensor<int, 2> mat2({5, 6, 7, 8}, 2, 2);
     tensor<int, 2> expected({6, 8, 10, 12}, 2, 2);

     auto result = mat1 + mat2;
     assert_matrix_equal(result, expected, "Addition Test");
 }

 void test_multiplication() {
     tensor<int, 2> mat1({1, 2, 3, 4, 5, 6}, 2, 3);
     tensor<int, 2> mat2({1, 1, 1, 1, 1, 1}, 3, 2);
     tensor<int, 2> expected({6, 6, 15, 15}, 2, 2);

     auto result = mat1.multiply(mat2);
     assert_matrix_equal(result, expected, "Multiplication Test");
 }

 bool floats_are_close(float a, float b, float tolerance = 0.01) {
     return std::fabs(a - b) < tolerance;
 }



void test_eigens()
{
    std::vector<double> data = {
        5,  10,  -5,
        2,  -14, 2,
        -4, -8,  6
   };
    squareMatrix<double> M(3, data);

    auto evals = M.eigenvalues();
    std::cout << "Eigenvalues:\n";
    for (auto &val : evals) {
        std::cout << val << "\n";
    }

    auto evecs = M.eigenvectors();
    std::cout << "Eigenvectors matrix (columns = eigenvectors):\n"
              << evecs.toString() << std::endl;
}

void test_inverse() {
    std::vector<double> data = {
            4, 7,
            2, 6
    };
    squareMatrix<double> A(2, data);
    auto Ainv = A.inverse();
    auto result = A.multiply(Ainv);

    // Generate identity matrix for comparison
    squareMatrix<double> expected = squareMatrix<double>::identity(2);

    bool pass = true;
    auto res_data = result.get_data();
    auto exp_data = expected.get_data();
    for (size_t i = 0; i < res_data.size(); ++i) {
        if (!floats_are_close(res_data[i], exp_data[i], 1e-6)) {
            std::cout << "FAIL [Inverse Test]: At index " << i << ", expected " << exp_data[i] << " but got " << res_data[i] << std::endl;
            pass = false;
        }
    }

    if (pass) {
        std::cout << "PASS [Inverse Test]" << std::endl;
    }
}

void test_lu_decomposition() {
    std::vector<double> data = {
            4, 3,
            6, 3
    };
    squareMatrix<double> A(2, data);

    auto [L, U] = A.luDecompose();
    auto recomposed = L.multiply(U);

    bool pass = true;
    auto orig = A.get_data();
    auto recomposed_data = recomposed.get_data();
    for (size_t i = 0; i < orig.size(); ++i) {
        if (!floats_are_close(orig[i], recomposed_data[i], 1e-6)) {
            std::cout << "FAIL [LU Decomposition Test]: At index " << i << ", expected " << orig[i] << " but got " << recomposed_data[i] << std::endl;
            pass = false;
        }
    }

    if (pass) {
        std::cout << "PASS [LU Decomposition Test]" << std::endl;
    }
}


void test_solve_equation() {
    std::vector<double> data = {
            2, 1,
            5, 7
    };
    squareMatrix<double> A(2, data);
    std::vector<double> b = {11, 13};

    std::vector<double> expected_x = {7.1111, -3.2222};

    auto x = A.solve(b);

    bool pass = true;
    for (size_t i = 0; i < x.size(); ++i) {
        if (!floats_are_close(x[i], expected_x[i], 0.01)) {
            std::cout << "FAIL [Solve Equation Test]: At index " << i
                      << ", expected " << expected_x[i] << " but got " << x[i] << std::endl;
            pass = false;
        }
    }

    if (pass) {
        std::cout << "PASS [Solve Equation Test]" << std::endl;
    }
}




int main() {
    test_addition();
    test_multiplication();
    test_eigens();
    test_inverse();
    test_lu_decomposition();
    test_solve_equation();
    std::cout << "All tests passed!" << std::endl;
    return 0;
}
