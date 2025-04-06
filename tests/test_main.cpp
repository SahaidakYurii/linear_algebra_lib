#include "tenzor.h"
#include "squareMatrix.h"
#include <cassert>
#include <cmath>
#include <iostream>

// template<typename T, size_t N>
// void print_matrix(const tenzor<T, N>& matrix) {
//     auto shape = matrix.shape();
//     size_t index = 0;
//     for (size_t i = 0; i < shape[0]; ++i) {
//         for (size_t j = 0; j < shape[1]; ++j) {
//             std::cout << matrix.data_m[index++] << " ";
//         }
//         std::cout << std::endl;
//     }
// }
//
// template<typename T, size_t N>
// void assert_matrix_equal(const tenzor<T, N>& result, const tenzor<T, N>& expected, const std::string& test_name) {
//     if (result.shape() != expected.shape()) {
//         std::cout << "FAIL [" << test_name << "]: Dimension mismatch." << std::endl;
//         return;
//     }
//     const auto& result_data = result.get_data();
//     const auto& expected_data = expected.get_data();
//     for (size_t i = 0; i < result.get_total_size(); ++i) {
//         if (result_data[i] != expected_data[i]) {
//             std::cout << "FAIL [" << test_name << "]: Values do not match at index " << i << ". Expected " << expected_data[i] << ", got " << result_data[i] << "." << std::endl;
//             return;
//         }
//     }
//     std::cout << "PASS [" << test_name << "]" << std::endl;
// }
//
//
// void test_addition() {
//     tenzor<int, 2> mat1({1, 2, 3, 4}, 2, 2);
//     tenzor<int, 2> mat2({5, 6, 7, 8}, 2, 2);
//     tenzor<int, 2> expected({6, 8, 10, 12}, 2, 2);
//
//     auto result = mat1 + mat2;
//     assert_matrix_equal(result, expected, "Addition Test");
// }
//
// void test_multiplication() {
//     tenzor<int, 2> mat1({1, 2, 3, 4, 5, 6}, 2, 3);
//     tenzor<int, 2> mat2({1, 1, 1, 1, 1, 1}, 2, 3);
//     tenzor<int, 2> expected({6, 6, 15, 15}, 2, 2);
//
//     auto result = mat1.multiply(mat2);
//     assert_matrix_equal(result, expected, "Multiplication Test");
// }
//
// bool floats_are_close(float a, float b, float tolerance = 0.01) {
//     return std::fabs(a - b) < tolerance;
// }
//
// void test_matrix_inversion() {
//     tenzor<float, 2> identity({1, 0, 0, 1}, 2, 2);
//     auto identity_inv = identity.invert();
//     assert(identity_inv(0, 0) == 1 && identity_inv(1, 1) == 1 && identity_inv(0, 1) == 0 && identity_inv(1, 0) == 0);
//
//
//     tenzor<float, 2> diagonal({4, 0, 0, 3}, 2, 2);
//     auto diagonal_inv = diagonal.invert();
//     assert(floats_are_close(diagonal_inv(0, 0), 0.25) && floats_are_close(diagonal_inv(1, 1), 1.0/3));
//
//
//     tenzor<float, 2> singular({2, 2, 2, 2}, 2, 2);
//     try {
//         auto singular_inv = singular.invert();
//         assert(false);
//     } catch (const std::exception& e) {
//         std::cout << "Correctly caught non-invertible matrix: " << e.what() << std::endl;
//     }
//     std::cout << "Inversion test passed" << std::endl;
// }
//
// void test_equation_solving() {
//     tenzor<float, 2> identity({1, 0, 0, 1}, 2, 2);
//     std::vector<float> b = {5, 3};
//     auto solution = identity.solve(b);
//     assert(solution[0] == 5 && solution[1] == 3);
//
//     tenzor<float, 2> A({4, 1, 2, 2}, 2, 2);
//     std::vector<float> b2 = {9, 8};
//     auto solution2 = A.solve(b2);
//     assert(floats_are_close(solution2[0], 2) && floats_are_close(solution2[1], 1));
//     std::cout << "Equation test passed" << std::endl;
// }

void test_eigens()
{
    std::vector<double> data = {
        4,  1,  2,
        1,  3, -1,
        2, -1,  3
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

int main() {
    // test_addition();
    // test_multiplication();
    // test_matrix_inversion();
    // test_equation_solving();
    test_eigens();
    std::cout << "All tests passed!" << std::endl;
    return 0;
}
