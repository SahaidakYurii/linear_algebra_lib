// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H
#include <sstream>

#include "tensor.h"

namespace linalg {
    template <typename T>
    class matrix : public tensor<T, 2> {
    protected:
        size_t rows_m, cols_m;

    public:
        matrix(const size_t rows, const size_t cols, vector<T> data = {})
            : tensor<T, 2>(data, rows, cols),
              rows_m(rows), cols_m(cols) {}

        explicit matrix(vector<vector<T>> vec2d)
            : tensor<T, 2>(
                vec2d.size(),
                vec2d.is_empty() ? 0 : std::max_element(vec2d.begin(), vec2d.end(),
                    [](const auto& a, const auto& b) {
                        return a.size() < b.size();
                    })->size()),
              rows_m(vec2d.size()),
              cols_m(vec2d.is_empty() ? 0 : std::max_element(vec2d.begin(), vec2d.end(),
                    [](const auto& a, const auto& b) {
                        return a.size() < b.size();
                    })->size())
        {
            for (size_t r = 0; r < rows_m; ++r) {
                std::copy(
                    vec2d[r].begin(),
                    vec2d[r].end(),
                    this->data_m.begin() + r * cols_m
                );

                if (vec2d[r].size() < cols_m) {
                    std::fill(
                        this->data_m.begin() + r * cols_m + vec2d[r].size(),
                        this->data_m.begin() + (r + 1) * cols_m,
                        T()
                    );
                }
            }
        }

        [[nodiscard]] size_t rows() const { return rows_m; }
        [[nodiscard]] size_t cols() const { return cols_m; }

        void reshape(const size_t rows, const size_t cols) {
            if (rows == rows_m && cols == cols_m) {
                return;
            }

            auto temp_data = vector<T>(rows * cols);

            for (size_t r = 0; r < (rows < rows_m ? rows : rows_m); r++) {
                for (size_t c = 0; c < (cols < cols_m ? cols : cols_m); c++) {
                    temp_data[r * cols + c] = this->data_m[r * cols_m + c];
                }
            }
            rows_m = rows; cols_m = cols;
            this->data_m = temp_data;
        }

        std::string toString() {
            std::stringstream ss;
            for (size_t r = 0; r < rows_m; r++) {
                for (size_t c = 0; c < cols_m; c++) {
                    ss << (*this)(r, c) << '\t';
                }
                ss << '\n';
            }

            return ss.str();
        }

        matrix<T> transpose() const{
            matrix<T> temp(cols_m, rows_m);
            for (size_t r = 0; r < rows_m; r++) {
                for (size_t c = 0; c < cols_m; c++) {
                    temp(c, r) = (*this)(r, c);
                }
            }

            return temp;
        }

        matrix<T>& operator*=(const matrix<T>& other) {
            if (this->cols_m != other.rows_m) {
                throw std::invalid_argument("Wrong dimensions of matrices");
            }

            size_t rows = this->rows_m, cols = other.cols_m;
            matrix<T> temp(rows, cols);

            for (size_t r = 0; r < rows; r++) {
                for (size_t c = 0; c < cols; c++) {
                    T sum{};
                    for (size_t i = 0; i < this->cols_m; i++) {
                        sum += (*this)(r, i) * other(i, c);
                    }
                    temp(r, c) = sum;
                }
            }
            *this = temp;
            return *this;
        }

        matrix<T>& operator*=(const T& scalar) {
            for (size_t i = 0; i < this->data_m.size(); ++i) {
                this->data_m[i] *= scalar;
            }
            return *this;
        }

        matrix<T> pinv() const;

        [[nodiscard]] size_t rank() const {
            matrix<T> temp(*this);
            const T EPS = static_cast<T>(1e-9);
            size_t rank = 0;
            size_t row = 0, col = 0;

            while (row < temp.rows_m && col < temp.cols_m) {
                // Find pivot
                size_t pivot = row;
                for (size_t i = row + 1; i < temp.rows_m; ++i) {
                    if (std::abs(temp(i, col)) > std::abs(temp(pivot, col))) {
                        pivot = i;
                    }
                }

                if (std::abs(temp(pivot, col)) < EPS) {
                    ++col;  // No pivot in this column
                    continue;
                }

                // Swap current row with pivot row
                for (size_t i = 0; i < temp.cols_m; ++i) {
                    std::swap(temp(row, i), temp(pivot, i));
                }

                // Eliminate below
                for (size_t i = row + 1; i < temp.rows_m; ++i) {
                    T factor = temp(i, col) / temp(row, col);
                    for (size_t j = col; j < temp.cols_m; ++j) {
                        temp(i, j) -= factor * temp(row, j);
                    }
                }

                ++rank;
                ++row;
                ++col;
            }

            return rank;
        }

        [[nodiscard]] vector<vector<T>> basis() const {
            matrix<T> temp(*this);
            const T EPS = static_cast<T>(1e-9);
            size_t m = temp.rows_m;
            size_t n = temp.cols_m;
            vector<size_t> pivot_columns;
            size_t row = 0;

            for (size_t col = 0; col < n && row < m; ++col) {
                size_t pivot = row;
                for (size_t i = row + 1; i < m; ++i) {
                    if (std::abs(temp(i, col)) > std::abs(temp(pivot, col))) {
                        pivot = i;
                    }
                }

                if (std::abs(temp(pivot, col)) < EPS) {
                    continue;
                }

                for (size_t j = 0; j < n; ++j) {
                    std::swap(temp(row, j), temp(pivot, j));
                }

                T pivot_val = temp(row, col);
                for (size_t j = col; j < n; ++j) {
                    temp(row, j) /= pivot_val;
                }

                for (size_t i = 0; i < m; ++i) {
                    if (i != row && std::abs(temp(i, col)) > EPS) {
                        T factor = temp(i, col);
                        for (size_t j = col; j < n; ++j) {
                            temp(i, j) -= factor * temp(row, j);
                        }
                    }
                }

                pivot_columns.push_back(col);
                ++row;
            }

            vector<vector<T>> basis_vectors;
            for (size_t col : pivot_columns) {
                vector<T> vec(m);
                for (size_t r = 0; r < m; ++r) {
                    vec[r] = (*this)(r, col);
                }
                basis_vectors.push_back(vec);
            }

            return basis_vectors;
        }

    }; // class matrix

    template <typename T>
    matrix<T> operator*(matrix<T> fst, const matrix<T>& snd) {
        fst *= snd;
        return fst;
    }

    template <typename T>
    matrix<T> operator*(matrix<T> fst, const T& snd) {
        fst *= snd;
        return fst;
    }

    template <typename T>
    matrix<T> operator*(const T& fst, matrix<T> snd) {
        snd *= fst;
        return snd;
    }

    template <typename T>
    matrix<T> operator+(const matrix<T>& lhs, const matrix<T>& rhs) {
        matrix<T> result(lhs);
        result += rhs;
        return result;
    }
} // namespace linalg

#include "squareMatrix.h"

namespace linalg {
    // Tell the compiler that there *will* be a class template squareMatrix<T>
    template<typename U>
    class squareMatrix;
}

namespace linalg{
    template <typename T>
    matrix<T> matrix<T>::pinv() const {
        auto AtA = static_cast<squareMatrix<T>>(this->transpose() * (*this));
        return AtA.inverse() * this->transpose();
    }
} // namespace linalg

#endif //LINALG_MATRIX_H
