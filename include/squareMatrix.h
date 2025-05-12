// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com


#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H

#include "matrix.h"

namespace linalg {
    template <typename T>
    class squareMatrix : public matrix<T> {
        size_t a;
    protected:
        squareMatrix<T> submatrix(size_t row, size_t col) {
            vector<T> sub_data{};
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

        vector<T> get_column(size_t colIndex) const {
            vector<T> col(this->a);
            for (size_t i = 0; i < this->a; i++) {
                col[i] = (*this)(i, colIndex);
            }
            return col;
        }

        void set_column(size_t colIndex, const vector<T>& colData)
        {
            for (size_t i = 0; i < this->a; i++) {
                (*this)(i, colIndex) = colData[i];
            }
        }

        T offDiagonalNorm() const {
            T sum = 0;
            for (size_t i = 0; i < a; i++) {
                for (size_t j = 0; j < a; j++) {
                    if (i != j) {
                        T val = (*this)(i,j);
                        sum += val * val;
                    }
                }
            }
            return std::sqrt(sum);
        }


    public:
        explicit squareMatrix(const size_t a, vector<T> data = {}):
            matrix<T>(a, a, data), a(a) {}

        explicit squareMatrix(const matrix<T> &other)
        : matrix<T>(other),  // or matrix<T>(other.data_m, other.rows_m, other.cols_m)
          a(other.rows_m)
        {
            if (other.rows_m != other.cols_m) {
                throw std::runtime_error("Not a square matrix!");
            }
        }

        explicit squareMatrix(matrix<T>&& other)
        : matrix<T>(std::move(other)), a(other.rows())
        {
            if (other.rows() != other.cols()) {
                throw std::runtime_error("Not a square matrix!");
            }
        }

        squareMatrix& operator=(const matrix<T>& other) {
            if (other.rows() != other.cols()) {
                throw std::runtime_error("Not a square matrix!");
            }

            if (this != &other) {
                matrix<T>::operator=(other);
                a = other.rows();
            }
            return *this;
        }

        squareMatrix& operator=(matrix<T>&& other) noexcept {
            if (other.rows() != other.cols()) {
                throw std::runtime_error("Not a square matrix!");
            }

            if (this != &other) {
                matrix<T>::operator=(std::move(other));
                a = other.rows();
            }

            return *this;
        }

        static squareMatrix<T> identity(size_t n) {
            vector<T> identity_data(n * n, T{0});
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
                    if ((*this)(r, c) != (*this)(c, r))
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

                    if ((*this)(r, c) != -(*this)(c, r))
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
                T val = (*this)(0, col);
                res += val * cof;
            }

            return res;
        }

        std::pair<squareMatrix<T>, squareMatrix<T>> qrDecompose() const {
            squareMatrix<T> Q(a);
            squareMatrix<T> R(a);

            for (size_t i = 0; i < a; i++) {
                for (size_t j = 0; j < a; j++) {
                    Q(i, j) = 0;
                    R(i, j) = 0;
                }
            }

            for (size_t k = 0; k < a; k++) {
                vector<T> vk = this->get_column(k);

                for (size_t j = 0; j < k; j++) {
                    vector<T> qj = Q.get_column(j);
                    T rjk = 0;
                    for (size_t i = 0; i < a; i++) {
                        rjk += qj[i] * vk[i];
                    }
                    R(j, k) = rjk;
                    for (size_t i = 0; i < a; i++) {
                        vk[i] -= rjk * qj[i];
                    }
                }

                T norm = vectorNorm(vk);
                if (norm < 1e-15) {
                    R(k, k) = 0;
                    for (size_t i = 0; i < a; i++) {
                        Q(i, k) = 0;
                    }
                } else {
                    R(k, k) = norm;
                    for (size_t i = 0; i < a; i++) {
                        Q(i, k) = vk[i] / norm;
                    }
                }
            }

            return {Q, R};
        }

        vector<T> eigenvalues(double tol = 1e-9, int maxIter = 1000) const {
            squareMatrix<T> A(*this);

            int iter = 0;
            while (iter < maxIter) {
                auto [Q, R] = A.qrDecompose();

                tensor<T, 2> multiplied = R.multiply(Q);
                A = squareMatrix<T>(A.a, multiplied.get_data());

                if (A.offDiagonalNorm() < tol) {
                    break;
                }
                iter++;
            }

            vector<T> eigvals(A.a);
            for (size_t i = 0; i < A.a; i++) {
                eigvals[i] = A(i, i);
            }
            return eigvals;
        }

        vector<T> nullSpace(double tol) const {
            squareMatrix<T> A = *this;
            size_t n = this->a;

            for (size_t i = 0; i < n; ++i) {
                size_t pivotRow = i;
                for (size_t j = i + 1; j < n; ++j) {
                    if (std::abs(A(j, i)) > std::abs(A(pivotRow, i))) {
                        pivotRow = j;
                    }
                }

                if (std::abs(A(pivotRow, i)) < tol) continue;

                if (pivotRow != i) {
                    for (size_t k = 0; k < n; ++k) std::swap(A(i, k), A(pivotRow, k));
                }

                for (size_t j = i + 1; j < n; ++j) {
                    T factor = A(j, i) / A(i, i);
                    for (size_t k = i; k < n; ++k) {
                        A(j, k) -= factor * A(i, k);
                    }
                }
            }

            vector<T> x(n, 0);
            x[n - 1] = 1;

            for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
                T sum = 0;
                for (size_t j = i + 1; j < n; ++j) {
                    sum += A(i, j) * x[j];
                }
                if (std::abs(A(i, i)) < tol) {
                    x[i] = 1;
                } else {
                    x[i] = -sum / A(i, i);
                }
            }

            return x;
        }

        squareMatrix<T> eigenvectorsViaNullspace(double tol = 1e-9) const {
            vector<T> eigenvals = this->eigenvalues();
            vector<vector<T>> vecs;

            for (const T& lambda : eigenvals) {
                squareMatrix<T> shifted = *this;
                for (size_t i = 0; i < this->a; ++i) {
                    shifted(i, i) -= lambda;
                }

                vector<T> v = shifted.nullSpace(tol);

                if (vectorNorm(v) > tol) {
                    for (auto& x : v) x /= vectorNorm(v);
                    vecs.push_back(v);
                } else {
                    vecs.push_back(vector<T>(this->a, 0));
                }
            }

            squareMatrix<T> result(this->a);
            for (size_t col = 0; col < vecs.size(); ++col) {
                for (size_t row = 0; row < this->a; ++row) {
                    result(row, col) = vecs[col][row];
                }
            }

            return result;
        }


        squareMatrix<T> inverse() {
            T det = this->determinant();
            if (std::abs(det) < 1e-12) {
                throw std::runtime_error("Matrix is singular or nearly singular — cannot invert.");
            }

            vector<T> cofactors_data(this->rows_m * this->cols_m);

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
                    U(i, j) = (*this)(i, j) - sum;
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
                        L(j, i) = ((*this)(j, i) - sum) / U(i, i);
                    }
                }
            }

            return std::make_pair(L, U);
        }

        vector<T> solve(const vector<T>& b) const {
            if (b.size() != this->a) {
                throw std::runtime_error("Size of b must match matrix dimensions.");
            }

            auto [L, U] = this->luDecompose();

            vector<T> y(this->a);
            for (size_t i = 0; i < this->a; ++i) {
                T sum = 0;
                for (size_t j = 0; j < i; ++j) {
                    sum += L(i, j) * y[j];
                }
                y[i] = b[i] - sum;
            }

            vector<T> x(this->a);
            for (int i = static_cast<int>(this->a) - 1; i >= 0; --i) {
                T sum = 0;
                for (size_t j = i + 1; j < this->a; ++j) {
                    sum += U(i, j) * x[j];
                }
                x[i] = (y[i] - sum) / U(i, i);
            }

            return x;
        }

        vector<vector<T>> eigenvectorsViaSolving(vector<T> eigenvalues, double tol = 1e-9) const {
            vector<vector<T>> eigenvectors;

            for (T lambda : eigenvalues) {
                // Form A - lambda * I
                squareMatrix<T> shifted = *this;
                for (size_t i = 0; i < this->a; ++i) {
                    shifted(i, i) -= lambda;
                }

                // We'll try to find x such that (A - λI)x = 0
                // Let’s feed a random b and normalize the result
                vector<T> b(this->a, 1.0); // or random vector, if you prefer
                vector<T> x;

                try {
                    x = shifted.solve(b);  // Will fail if matrix is singular
                } catch (...) {
                    // Matrix likely singular => try alternate approach
                    x = vector<T>(this->a, 0);
                }

                // Normalize eigenvector
                T norm = vectorNorm(x);
                if (norm > tol) {
                    for (T& val : x) val /= norm;
                    eigenvectors.push_back(x);
                }
            }

            return eigenvectors;
        }
    };

    template <typename T>
    squareMatrix<T> operator+(matrix<T> lhs, const matrix<T>& rhs) {
        lhs += rhs;
        return lhs;
    }
} // namespace linalg

#endif //SQUAREMATRIX_H
