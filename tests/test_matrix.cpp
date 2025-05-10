#include <gtest/gtest.h>
#include "linalg.h"
#include <cmath>

TEST(MatrixTests, EigenvaluesAndVectors) {
    linalg::vector<double> data = {
        5, 10, -5,
        2, -14, 2,
        -4, -8, 6
    };
    linalg::squareMatrix<double> M(3, data);
    auto evals = M.eigenvalues();

    for (const auto& val : evals) {
        EXPECT_TRUE(std::isfinite(val));
    }

    auto evecs = M.eigenvectorsViaNullspace();
    EXPECT_EQ(evecs.rows(), 3);
    EXPECT_EQ(evecs.cols(), evals.size());
}

TEST(MatrixTests, Inverse) {
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

TEST(MatrixTests, LUDecomposition) {
    linalg::squareMatrix<double> A(2, {4, 3, 6, 3});
    auto [L, U] = A.luDecompose();
    auto recomposed = L.multiply(U);

    const auto& orig = A.get_data();
    const auto& recomposed_data = recomposed.get_data();
    for (size_t i = 0; i < orig.size(); ++i) {
        EXPECT_NEAR(orig[i], recomposed_data[i], 1e-6);
    }
}

TEST(MatrixTests, SolveEquation) {
    linalg::squareMatrix<double> A(2, {2, 1, 5, 7});
    linalg::vector<double> b = {11, 13};
    linalg::vector<double> expected_x = {7.1111, -3.2222};

    auto x = A.solve(b);

    ASSERT_EQ(x.size(), expected_x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(x[i], expected_x[i], 0.01);
    }
}

TEST(MatrixTests, PseudoInverse) {
    linalg::matrix<double> A(3, 2, {1, 2, 3, 4, 5, 6});
    auto A_p = A.pinv();

    // A * A.pinv() * A == A
    auto AA_pA = A * A_p * A;
    for (size_t r = 0; r < A.rows(); ++r) {
        for (size_t c = 0; c < A.cols(); ++c) {
            EXPECT_NEAR(A(r, c), AA_pA(r, c), 0.001);
        }
    }

    // A.pinv() * A * A.pinv() == A.pinv()
    auto A_pAA_p = A_p * A * A_p;
    for (size_t r = 0; r < A.rows(); ++r) {
        for (size_t c = 0; c < A.cols(); ++c) {
            EXPECT_NEAR(A_p(r, c), A_pAA_p(r, c), 0.001);
        }
    }

     // Check if A.pinv() * A is symmetrical
     auto ApA = A_p * A;
     ASSERT_EQ(ApA.rows(), ApA.cols()) << "A.pinv() * A is not square!";
     for (size_t i = 0; i < ApA.rows(); ++i) {
         for (size_t j = 0; j < ApA.cols(); ++j) {
             EXPECT_NEAR(ApA(i, j), ApA(j, i), 0.1);
         }
     }

     // Check if A * A.pinv() is symmetrical
     auto AAp = A * A_p;
     ASSERT_EQ(AAp.rows(), AAp.cols()) << "A * A.pinv() is not square!";
     for (size_t i = 0; i < AAp.rows(); ++i) {
         for (size_t j = 0; j < AAp.cols(); ++j) {
             EXPECT_NEAR(AAp(i, j), AAp(j, i), 0.1);
         }
     }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

