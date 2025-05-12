#ifndef LINALG_SPARSE_MATRIX_H
#define LINALG_SPARSE_MATRIX_H

#include <unordered_map>
#include <tuple>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <unordered_map>
#include <tuple>
#include <functional>
#include "vector.h"

namespace std {
    template <>
    struct hash<std::tuple<size_t, size_t>> {
        size_t operator()(const std::tuple<size_t, size_t>& key) const noexcept {
            auto [row, col] = key;
            return std::hash<size_t>()(row) ^ (std::hash<size_t>()(col) << 1);
        }
    };
}

namespace linalg {
    using coords = std::tuple<size_t, size_t>;

    template<typename T>
    class sparseMatrix {
    private:
        size_t rows_m, cols_m;
        std::unordered_map<coords, T> data_m;

    public:
        sparseMatrix(const size_t rows, const size_t cols, vector<T> data = {})
            : rows_m(rows), cols_m(cols), data_m{} {
            if (data.is_empty())
                return;

            if (data.size() != rows_m * cols_m) {
                throw std::invalid_argument("data must be of size rows * cols");
            }

            for (size_t i = 0; i < data.size(); i++) {
                if (data[i] != T())
                    data_m[{i / cols_m, i % cols_m}] = data[i];
            }
        }

        [[nodiscard]] size_t rows() const { return rows_m; }
        [[nodiscard]] size_t cols() const { return cols_m; }

        void reshape(const size_t rows, const size_t cols) {
            if (rows == rows_m && cols == cols_m) {
                return;
            }

            auto temp_data = std::unordered_map<coords, T>();

            for (const auto& [key, value] : data_m) {
                auto [row, col] = key;
                if (row < rows && col < cols) {
                    temp_data[{row, col}] = value;
                }
            }
            data_m = std::move(temp_data);
            rows_m = rows; cols_m = cols;
        }

        [[nodiscard]] std::string toString() const {
            std::stringstream ss;
            for (size_t r = 0; r < rows_m; r++) {
                for (size_t c = 0; c < cols_m; c++) {
                    if (data_m.contains({r, c}))
                        ss << data_m.at({r, c}) << '\t';
                    else
                        ss << T() << '\t';
                }
            }

            return ss.str();
        }

        sparseMatrix<T> transpose() const {
            sparseMatrix<T> temp(cols_m, rows_m);
            for (const auto& [key, value] : data_m) {
                auto [r, c] = key;
                temp(c, r) = value;
            }
            return  temp;
        }

        T& operator()(const size_t row, const size_t col) {
            return data_m[{row, col}];
        }

        T operator()(const size_t row, const size_t col) const {
            auto it = data_m.find({row, col});
            return (it != data_m.end()) ? it->second : T();
        }

        sparseMatrix<T> operator+=(const sparseMatrix<T>& other) const {
            for (const auto& [key, value] : other.data_m) {
                auto [r, c] = key;
                data_m[{r, c}] += value;
            }
            return *this;
        }

        sparseMatrix<T>& operator*=(const sparseMatrix<T>& other) {
            if (this->cols_m != other.rows_m) {
                throw std::invalid_argument("Wrong dimensions of matrices");
            }

            size_t rows = this->rows_m, cols = other.cols_m;
            sparseMatrix<T> temp(rows, cols);

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

        sparseMatrix<T>& operator*=(const T& scalar) {
            for (const auto& [key, value] : data_m) {
                this->data_m[key] = scalar * value;
            }

            return *this;
        }

    }; // class sparseMatrix

    template <typename T>
    sparseMatrix<T> operator+(sparseMatrix<T> fst, const sparseMatrix<T>& snd) {
        fst += snd;
        return fst;
    }

    template <typename T>
    sparseMatrix<T> operator*(sparseMatrix<T> fst, const sparseMatrix<T>& snd) {
        fst *= snd;
        return fst;
    }

    template <typename T>
    sparseMatrix<T> operator*(sparseMatrix<T> fst, const T& snd) {
        fst *= snd;
        return fst;
    }

    template <typename T>
    sparseMatrix<T> operator*(const T& fst, sparseMatrix<T> snd) {
        snd *= fst;
        return snd;
    }
} // namespace linalg

#endif // LINALG_SPARSE_MATRIX_H