#include <gtest/gtest.h>
#include "linalg.h"
#include <cmath>

bool floats_are_close(double a, double b, double tolerance = 0.01) {
    return std::fabs(a - b) < tolerance;
}

template<typename T, size_t N>
void assert_tensor_equal(const linalg::tensor<T, N>& result, const linalg::tensor<T, N>& expected) {
    ASSERT_EQ(result.shape(), expected.shape()) << "Shape mismatch.";

    const auto& result_data = result.get_data();
    const auto& expected_data = expected.get_data();
    ASSERT_EQ(result_data.size(), expected_data.size());

    for (size_t i = 0; i < result_data.size(); ++i) {
        EXPECT_EQ(result_data[i], expected_data[i]) << "Mismatch at index " << i;
    }
}

TEST(TensorTests, Addition) {
    linalg::tensor<int, 2> mat1({1, 2, 3, 4}, 2, 2);
    linalg::tensor<int, 2> mat2({5, 6, 7, 8}, 2, 2);
    linalg::tensor<int, 2> expected({6, 8, 10, 12}, 2, 2);

    auto result = mat1 + mat2;
    assert_tensor_equal(result, expected);
}

TEST(TensorTests, Multiplication) {
    linalg::tensor<int, 2> mat1({1, 2, 3, 4, 5, 6}, 2, 3);
    linalg::tensor<int, 2> mat2({1, 1, 1, 1, 1, 1}, 3, 2);
    linalg::tensor<int, 2> expected({6, 6, 15, 15}, 2, 2);

    auto result = mat1.multiply(mat2);
    assert_tensor_equal(result, expected);
}
