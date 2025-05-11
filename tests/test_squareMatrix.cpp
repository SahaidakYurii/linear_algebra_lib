#include <gtest/gtest.h>
#include "linalg.h"

TEST(SquareMatrixTests, Inverse) {
    linalg::squareMatrix<double> A(2, {4, 7, 2, 6});
    auto Ainv = A.inverse();
    auto result = A.multiply(Ainv);
    auto expected = linalg::squareMatrix<double>::identity(2);

    const auto& res_data = result.get_data();
    const auto& exp_data = expected.get_data();

    for (size_t i = 0; i < res_data.size(); ++i) {
        EXPECT_NEAR(res_data[i], exp_data[i], 1e-6);
    }
}

TEST(SquareMatrixTests, LUDecomposition) {
    linalg::squareMatrix<double> A(2, {4, 3, 6, 3});
    auto [L, U] = A.luDecompose();
    auto recomposed = L.multiply(U);

    const auto& orig = A.get_data();
    const auto& recomposed_data = recomposed.get_data();
    for (size_t i = 0; i < orig.size(); ++i) {
        EXPECT_NEAR(orig[i], recomposed_data[i], 1e-6);
    }
}

TEST(SquareMatrixTests, SolveEquation) {
    linalg::squareMatrix<double> A(2, {2, 1, 5, 7});
    linalg::vector<double> b = {11, 13};
    linalg::vector<double> expected_x = {7.1111, -3.2222};

    auto x = A.solve(b);

    ASSERT_EQ(x.size(), expected_x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(x[i], expected_x[i], 0.01);
    }
}