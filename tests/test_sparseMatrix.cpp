#include <gtest/gtest.h>
#include "sparseMatrix.h"

using namespace linalg;

TEST(SparseMatrixTest, ConstructionAndAccess) {
    sparseMatrix<double> m(3, 3,{
        0, 1, 0,
        0, 0, 2,
        0, 0, 3});
    EXPECT_DOUBLE_EQ(m(0, 1), 1);
    EXPECT_DOUBLE_EQ(m(1, 2), 2);
    EXPECT_DOUBLE_EQ(m(2, 2), 3);
    EXPECT_DOUBLE_EQ(m(0, 0), 0);
}

TEST(SparseMatrixTest, Reshape) {
    sparseMatrix<int> m(3, 3);
    m(0, 1) = 5;
    m(2, 2) = 7;
    m.reshape(2, 2);
    EXPECT_EQ(m.rows(), 2);
    EXPECT_EQ(m.cols(), 2);
    EXPECT_EQ(m(0, 1), 5);
    EXPECT_EQ(m(1, 1), 0);
}

TEST(SparseMatrixTest, Transpose) {
    sparseMatrix<int> m(2, 3);
    m(0, 1) = 10;
    m(1, 2) = 20;
    auto t = m.transpose();
    EXPECT_EQ(t.rows(), 3);
    EXPECT_EQ(t.cols(), 2);
    EXPECT_EQ(t(1, 0), 10);
    EXPECT_EQ(t(2, 1), 20);
    EXPECT_EQ(t(0, 0), 0);
}

TEST(SparseMatrixTest, ScalarMultiplication) {
    sparseMatrix<float> m(2, 2);
    m(0, 0) = 1.5f;
    m(1, 1) = 2.5f;
    m *= 2.0f;
    EXPECT_FLOAT_EQ(m(0, 0), 3.0f);
    EXPECT_FLOAT_EQ(m(1, 1), 5.0f);
    EXPECT_FLOAT_EQ(m(0, 1), 0.0f);
}

TEST(SparseMatrixTest, MatrixMultiplication) {
    sparseMatrix<int> a(2, 3);
    a(0, 1) = 1;
    a(1, 2) = 2;

    sparseMatrix<int> b(3, 2);
    b(1, 0) = 3;
    b(2, 1) = 4;

    auto c = a * b;
    EXPECT_EQ(c.rows(), 2);
    EXPECT_EQ(c.cols(), 2);
    EXPECT_EQ(c(0, 0), 3);
    EXPECT_EQ(c(1, 1), 8);
    EXPECT_EQ(c(0, 1), 0);
}

TEST(SparseMatrixTest, MatrixMultiplicationDataRace) {
    constexpr size_t N = 100;
    sparseMatrix<int> a(1, N);
    sparseMatrix<int> b(N, N);

    for (size_t i = 0; i < N; ++i)
        a(0, i) = 1;

    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j)
            b(i, j) = 1;

    auto c = a * b;

    for (size_t i = 1; i < N; ++i) {
        EXPECT_EQ(c(0, i), 100);
    }
}