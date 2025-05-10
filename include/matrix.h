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

        [[nodiscard]] size_t rows() const { return rows_m; }
        [[nodiscard]] size_t cols() const { return cols_m; }

        void reshape(const size_t rows, const size_t cols) {
            if (rows == rows_m || cols == cols_m) {
                return;
            }

            auto temp_data = vector<T>(rows * cols);

            // copies the common part, if 3x4 reshaped to 2x6, the 2x4 part will remain same
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
                    ss << this->data_m[r * cols_m + c] << '\t';
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
    };

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
namespace linalg{
    template <typename T>
    matrix<T> matrix<T>::pinv() const {
        auto AtA = static_cast<squareMatrix<T>>(this->transpose() * (*this));
        return AtA.inverse() * this->transpose();
    }
} // namespace linalg

#endif //LINALG_MATRIX_H
