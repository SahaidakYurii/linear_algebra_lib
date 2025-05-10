#include <gtest/gtest.h>
#include "linalg.h"
#include <cmath>

bool floats_are_close(double a, double b, double tolerance = 0.01) {
    return std::fabs(a - b) < tolerance;
}

template<typename T, size_t N>
void assert_tensor_equal(const tensor<T, N>& result, const tensor<T, N>& expected) {
    ASSERT_EQ(result.shape(), expected.shape()) << "Shape mismatch.";

    const auto& result_data = result.get_data();
    const auto& expected_data = expected.get_data();
    ASSERT_EQ(result_data.size(), expected_data.size());

    for (size_t i = 0; i < result_data.size(); ++i) {
        EXPECT_EQ(result_data[i], expected_data[i]) << "Mismatch at index " << i;
    }
}

TEST(TensorTests, Addition) {
    tensor<int, 2> mat1({1, 2, 3, 4}, 2, 2);
    tensor<int, 2> mat2({5, 6, 7, 8}, 2, 2);
    tensor<int, 2> expected({6, 8, 10, 12}, 2, 2);

    auto result = mat1 + mat2;
    assert_tensor_equal(result, expected);
}

TEST(TensorTests, Multiplication) {
    tensor<int, 2> mat1({1, 2, 3, 4, 5, 6}, 2, 3);
    tensor<int, 2> mat2({1, 1, 1, 1, 1, 1}, 3, 2);
    tensor<int, 2> expected({6, 6, 15, 15}, 2, 2);

    auto result = mat1.multiply(mat2);
    assert_tensor_equal(result, expected);
}

TEST(MatrixTests, EigenvaluesAndVectors) {
    std::vector<double> data = {
        5, 10, -5,
        2, -14, 2,
        -4, -8, 6
    };
    squareMatrix<double> M(3, data);
    auto evals = M.eigenvalues();

    for (const auto& val : evals) {
        EXPECT_TRUE(std::isfinite(val));
    }

    auto evecs = M.eigenvectorsViaNullspace();
    EXPECT_EQ(evecs.rows(), 3);
    EXPECT_EQ(evecs.cols(), evals.size());
}

TEST(MatrixTests, Inverse) {
    squareMatrix<double> A(2, {4, 7, 2, 6});
    auto Ainv = A.inverse();
    auto result = A.multiply(Ainv);
    auto expected = squareMatrix<double>::identity(2);

    const auto& res_data = result.get_data();
    const auto& exp_data = expected.get_data();

    for (size_t i = 0; i < res_data.size(); ++i) {
        EXPECT_NEAR(res_data[i], exp_data[i], 1e-6);
    }
}

TEST(MatrixTests, LUDecomposition) {
    squareMatrix<double> A(2, {4, 3, 6, 3});
    auto [L, U] = A.luDecompose();
    auto recomposed = L.multiply(U);

    const auto& orig = A.get_data();
    const auto& recomposed_data = recomposed.get_data();
    for (size_t i = 0; i < orig.size(); ++i) {
        EXPECT_NEAR(orig[i], recomposed_data[i], 1e-6);
    }
}

TEST(MatrixTests, SolveEquation) {
    squareMatrix<double> A(2, {2, 1, 5, 7});
    std::vector<double> b = {11, 13};
    std::vector<double> expected_x = {7.1111, -3.2222};

    auto x = A.solve(b);

    ASSERT_EQ(x.size(), expected_x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(x[i], expected_x[i], 0.01);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

