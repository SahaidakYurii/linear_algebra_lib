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

    auto AAp = A * A_p;
    ASSERT_EQ(AAp.rows(), AAp.cols()) << "A * A.pinv() is not square!";
    for (size_t i = 0; i < AAp.rows(); ++i) {
        for (size_t j = 0; j < AAp.cols(); ++j) {
            EXPECT_NEAR(AAp(i, j), AAp(j, i), 0.1);
        }
    }
}

TEST(MatrixTests, RankCalculation) {
    using linalg::matrix;

    matrix<double> A(3, 3, {
        1, 2, 3,
        0, 1, 4,
        0, 0, 5
    });
    EXPECT_EQ(A.rank(), 3);

    matrix<double> B(3, 3, {
        1, 2, 3,
        2, 4, 6,
        3, 6, 9
    });
    EXPECT_EQ(B.rank(), 1);

    matrix<double> C(3, 3, {
        1, 2, 3,
        4, 5, 6,
        0, 0, 0
    });
    EXPECT_EQ(C.rank(), 2);

    matrix<double> D(4, 4, {
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    });
    EXPECT_EQ(D.rank(), 4);

    matrix<double> E(4, 2, {
        1, 0,
        0, 1,
        1, 1,
        2, 1
    });
    EXPECT_EQ(E.rank(), 2);

    matrix<double> F(2, 4, {
        1, 2, 3, 4,
        2, 4, 6, 8
    });
    EXPECT_EQ(F.rank(), 1);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

