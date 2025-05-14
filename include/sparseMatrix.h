#ifndef LINALG_SPARSE_MATRIX_H
#define LINALG_SPARSE_MATRIX_H

#include <unordered_map>
#include <tuple>
#include <stdexcept>
#include <cmath>
#include <sstream>
#include <functional>
#include "vector.h"
#include "linalg.h"

#if LINALG_USE_THREADS
  #include <tbb/concurrent_unordered_map.h>
  using tbb::concurrent_unordered_map;
#else
    #include <unordered_map>
#endif

template <>
struct std::hash<std::tuple<size_t, size_t>> {
    size_t operator()(const std::tuple<size_t, size_t>& key) const noexcept {
        auto [row, col] = key;
        size_t seed = std::hash<size_t>{}(row);
        seed ^= std::hash<size_t>{}(col)
                + 0x9e3779b97f4a7c15ULL
                + (seed << 6)
                + (seed >> 2);
        return seed;
    }
};

namespace linalg {
    using coords = std::tuple<size_t, size_t>;

    template<typename T>
    class sparseMatrix {
    private:
        size_t rows_m, cols_m;
#if LINALG_USE_THREADS
        concurrent_unordered_map<coords, T> data_m;
#else
        std::unordered_map<coords, T> data_m;
#endif

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

            auto temp_data = decltype(data_m){};

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

        T& operator()(size_t r, size_t c) {
            return data_m[{r,c}];
        }
        T operator()(size_t r, size_t c) const {
            auto it = data_m.find({r,c});
            return it != data_m.end() ? it->second : T();
        }

        sparseMatrix<T> operator+=(const sparseMatrix<T>& other) {
            for (const auto& [key, value] : other.data_m) {
                auto [r, c] = key;
                data_m[{r, c}] += value;
            }
            return *this;
        }

        sparseMatrix<T>& operator*=(const sparseMatrix<T>& other) {
            if (cols_m != other.rows_m)
                throw std::invalid_argument("Wrong dimensions");

            size_t R = rows_m, C = other.cols_m;
#if LINALG_USE_THREADS
            concurrent_unordered_map<coords, T> temp;
            detail::ThreadPool& pool = detail::ThreadPool::instance();
            std::vector<std::future<void>> tasks;

            for (auto const& [key, v] : data_m) {
                auto [r,c] = key;
                tasks.emplace_back(pool.enqueue([&,r,c,v]{
                    for (size_t j = 0; j < C; ++j) {
                        T ov = other(c,j);
                        if (ov != T()) {
                            temp[{r,j}] += v * ov;
                        }
                    }
                }));
            }
            for (auto &f : tasks) f.get();

            sparseMatrix<T> result(R, C);
            for (auto const& kv : temp)
                result.data_m.insert(kv);
            *this = std::move(result);

#else
            sparseMatrix<T> result(R, C);
            for (auto const& [key, v] : data_m) {
                auto [r,c] = key;
                for (size_t j = 0; j < C; ++j) {
                    T ov = other(c,j);
                    if (ov != T())
                        result(r,j) += v * ov;
                }
            }
            *this = std::move(result);
#endif

            return *this;
        }

        sparseMatrix<T>& operator*=(const T& scalar) {
            for (const auto& [key, value] : data_m) {
                this->data_m[key] = scalar * value;
            }

            return *this;
        }

    };

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