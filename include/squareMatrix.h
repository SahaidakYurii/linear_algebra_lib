//
// Created by yurii-sahaidak on 4/5/25.
//

#include "matrix.h"

#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H

template <typename T>
class squareMatrix : public matrix<T> {
    size_t a;
protected:
    squareMatrix<T> submatrix(size_t row, size_t col) {
        std::vector<T> sub_data{};
        sub_data.reserve((this->rows_m - 1) * (this->cols_m - 1));
        for (size_t r = 0; r < this->rows_m; r++) {
            if (r == row) continue;
            for (size_t c = 0; c < this->cols_m; c++) {
                if (c == col) continue;

                sub_data.push_back(this->data_m[r * this->cols_m + c]);
            }
        }

        return squareMatrix<T>{a-1, sub_data};
    }

    T vectorNorm(const std::vector<T>& v) const {
        T sum = 0;
        for (auto &val : v) sum += val * val;
        return std::sqrt(sum);
    }

    std::vector<T> get_column(const squareMatrix<T>& mat, size_t colIndex) const {
        std::vector<T> col(mat.a);
        for (size_t i = 0; i < mat.a; i++) {
            col[i] = mat(i, colIndex);
        }
        return col;
    }

    static void set_column(squareMatrix<T>& mat, size_t colIndex, const std::vector<T>& colData)
    {
        for (size_t i = 0; i < mat.a; i++) {
            mat.data_m[mat.a + i] = colData[i];
        }
    }

    std::pair<squareMatrix<T>, squareMatrix<T>> qrDecompose() const {
        squareMatrix<T> Q(a);
        squareMatrix<T> R(a);
        for (size_t i = 0; i < a; i++) {
            for (size_t j = 0; j < a; j++) {
                Q(i, j) = (i == j) ? 1 : 0;
               R(i, j) = 0;

                // Q.data_m[i * a + j] = (i == j) ? 1 : 0;
                // R.data_m[i * a + j] = 0;
            }
        }

        for (size_t k = 0; k < a; k++) {
            std::vector<T> col = get_column(*this, k);

            for (size_t j = 0; j < k; j++) {
                std::vector<T> qj = get_column(Q, j);
                T dot = 0;
                for (size_t idx = 0; idx < a; idx++) {
                    dot += qj[idx] * col[idx];
                }
                R(j, k) = dot;
                for (size_t idx = 0; idx < a; idx++) {
                    col[idx] -= dot * qj[idx];
                }
            }

            T norm_col = vectorNorm(col);
            if (norm_col < 1e-15) {
                R(k, k) = 0;
                for (auto &val : col) val = 0;
            } else {
                R(k, k) = norm_col;
                for (size_t idx = 0; idx < a; idx++) {
                    col[idx] /= norm_col;
                }
            }

            set_column(Q, k, col);
        }

        return std::make_pair(Q, R);
    }

    T offDiagonalNorm() const {
        T sum = 0;
        for (size_t i = 0; i < a; i++) {
            for (size_t j = 0; j < a; j++) {
                if (i != j) {
                    T val = this->operator()(i,j);
                    sum += val * val;
                }
            }
        }
        return std::sqrt(sum);
    }


public:
    explicit squareMatrix(const size_t a, std::vector<T> data = {}):
        matrix<T>(a, a, data), a(a) {}

    explicit squareMatrix(const matrix<T> &other)
    : matrix<T>(other),  // or matrix<T>(other.data_m, other.rows_m, other.cols_m)
      a(other.rows_m)
    {
        if (other.rows_m != other.cols_m) {
            throw std::runtime_error("Not a square matrix!");
        }
    }




    static squareMatrix<T> identity(size_t n) {
        std::vector<T> identity_data(n * n, T{0});
        for (size_t i = 0; i < n; ++i) {
            identity_data[i * n + i] = T{1};
        }
        return squareMatrix<T>(n, identity_data);
    }

    bool isSingular() {
        return determinant() == 0;
    }

    bool isSymmetrical() {
        for (size_t r = 0; r < a; r++) {
            for (size_t c = 0; c < a; c++) {
                if (this->operator()(r, c) != this->operator()(c, r))
                    return false;
            }
        }
        return true;
    }

    bool isSkewSymmetrical() {
        for (size_t r = 0; r < a; r++) {
            for (size_t c = 0; c < a; c++) {
                if (r == c)
                    continue;
                
                if (this->operator()(r, c) != -this->operator()(c, r))
                    return false;
            }
        }
        return true;
    }

    T minor(size_t row, size_t col) {
        squareMatrix<T> sub = submatrix(row, col);
        return sub.determinant();
    }

    T cofactor(size_t row, size_t col) {
        return (((row + col) % 2 == 0) ? 1 : -1) * minor(row, col);
    }

    T determinant() {
        if (this->rows_m == 1)
            return this->data_m[0];

        if (this->rows_m == 2)
            return this->data_m[0] * this->data_m[3] - this->data_m[1] * this->data_m[2];

        T res = 0;
        for (size_t col = 0; col < this->cols_m; ++col) {
            T cof = cofactor(0, col);
            T val = this->operator()(0, col);
            res += val * cof;
        }

        return res;
    }

    std::vector<T> eigenvalues(double tol = 1e-9, int maxIter = 1000) const {
        squareMatrix<T> A(*this);

        int iter = 0;
        while (iter < maxIter) {
            auto [Q, R] = A.qrDecompose();

            tenzor<T, 2> multiplied = R.multiply(Q);
            A = squareMatrix<T>(A.a, multiplied.get_data());

            if (A.offDiagonalNorm() < tol) {
                break;
            }
            iter++;
        }

        std::vector<T> eigvals(A.a);
        for (size_t i = 0; i < A.a; i++) {
            eigvals[i] = A(i, i);
        }
        return eigvals;
    }

    squareMatrix<T> eigenvectors(double tol = 1e-9, int maxIter = 1000) const {

        squareMatrix<T> A(*this);
        squareMatrix<T> Qaccum(A.a);
        for (size_t i = 0; i < A.a; i++) {
            for (size_t j = 0; j < A.a; j++) {
                Qaccum(i, j) = (i == j) ? 1 : 0;
            }
        }

        int iter = 0;
        while (iter < maxIter) {
            auto [Q, R] = A.qrDecompose();
            tenzor<T, 2> multiplied = R.multiply(Q);
            A = squareMatrix<T>(A.a, multiplied.get_data());
            tenzor<T, 2> qaccumMult = Qaccum.multiply(Q);
            Qaccum = squareMatrix<T>(A.a, qaccumMult.get_data());
            if (A.offDiagonalNorm() < tol) {
                break;
            }
            iter++;
        }
        return Qaccum;
    }

    squareMatrix<T> inverse() {
        T det = this->determinant();
        if (std::abs(det) < 1e-12) {
            throw std::runtime_error("Matrix is singular or nearly singular — cannot invert.");
        }

        std::vector<T> cofactors_data(this->rows_m * this->cols_m);

        for (size_t r = 0; r < this->rows_m; r++) {
            for (size_t c = 0; c < this->cols_m; c++) {
                cofactors_data[r * this->cols_m + c] = this->cofactor(r, c);
            }
        }

        squareMatrix<T> cofactorMat(this->rows_m, cofactors_data);

        squareMatrix<T> adjugate(this->rows_m);
        for (size_t i = 0; i < this->rows_m; ++i) {
            for (size_t j = 0; j < this->cols_m; ++j) {
                adjugate(i, j) = cofactorMat(j, i);
            }
        }

        for (auto& val : adjugate.data_m) {
            val /= det;
        }

        return adjugate;
    }

    std::pair<squareMatrix<T>, squareMatrix<T>> luDecompose() const {
        if (this->a != this->rows_m || this->a != this->cols_m) {
            throw std::runtime_error("LU decomposition only valid for square matrices");
        }

        squareMatrix<T> L(this->a);
        squareMatrix<T> U(this->a);

        for (size_t i = 0; i < this->a; i++) {
            // Fill U
            for (size_t j = i; j < this->a; j++) {
                T sum = 0;
                for (size_t k = 0; k < i; k++) {
                    sum += L(i, k) * U(k, j);
                }
                U(i, j) = this->operator()(i, j) - sum;
            }

            // Fill L
            for (size_t j = i; j < this->a; j++) {
                if (i == j) {
                    L(i, i) = 1;
                } else {
                    T sum = 0;
                    for (size_t k = 0; k < i; k++) {
                        sum += L(j, k) * U(k, i);
                    }
                    L(j, i) = (this->operator()(j, i) - sum) / U(i, i);
                }
            }
        }

        return std::make_pair(L, U);
    }

    std::vector<T> solve(const std::vector<T>& b) const {
        if (b.size() != this->a) {
            throw std::runtime_error("Size of b must match matrix dimensions.");
        }

        auto [L, U] = this->luDecompose();

        std::vector<T> y(this->a);
        for (size_t i = 0; i < this->a; ++i) {
            T sum = 0;
            for (size_t j = 0; j < i; ++j) {
                sum += L(i, j) * y[j];
            }
            y[i] = b[i] - sum;
        }

        std::vector<T> x(this->a);
        for (int i = static_cast<int>(this->a) - 1; i >= 0; --i) {
            T sum = 0;
            for (size_t j = i + 1; j < this->a; ++j) {
                sum += U(i, j) * x[j];
            }
            x[i] = (y[i] - sum) / U(i, i);
        }

        return x;
    }




};

#endif //SQUAREMATRIX_H
