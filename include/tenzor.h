// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#ifndef LINALG_NDMATRIX_H
#define LINALG_NDMATRIX_H

#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>
#include <initializer_list>
#include <cmath>

template <typename T, size_t N>
class tenzor {
protected:
    std::vector<T> data_m;
    std::array<size_t, N> dims_m;

    size_t total_size() const {
        size_t total = 1;
        for (size_t i = 0; i < N; ++i) {
            total *= dims_m[i];
        }
        return total;
    }

    template <typename... Indices>
    size_t calculate_index(size_t iter, size_t idx, Indices... indices) const {
        size_t stride = 1;
        for (size_t i = iter + 1; i < N; ++i) {
            stride *= dims_m[i];
        }
        return idx * stride + calculate_index(iter + 1, indices...);

    }
    size_t calculate_index(size_t iter) const {
        if (iter != N) {
            throw std::invalid_argument("Number of arguments does not match dimentionality of the tenzor");
        }
        return 0;
    }

public:
    template <typename... Indices>
    explicit tenzor(Indices ... indices) : data_m() {
        static_assert(sizeof...(Indices) == N, "Number of indices must match the dimension N");
        static_assert((std::is_convertible_v<Indices, size_t> && ...), "All indices must be convertible to size_t");

        dims_m = { static_cast<size_t>(indices)... };
        data_m.resize(total_size());
    }

    template <typename... Indices>
    explicit tenzor(const std::vector<T>& data, Indices... indices)
            : tenzor(indices...){
        data_m = data;
        data_m.resize(total_size());

    }

    tenzor<T, N>& reshape(std::initializer_list<size_t> new_dims) {
        if (new_dims.size() != N) {
            throw std::invalid_argument("Dimension size does not match template parameter N.");
        }

        size_t new_total_size = 1;
        for (auto dim : new_dims) {
            new_total_size *= dim;
        }

        if (new_total_size != total_size()) {
            throw std::invalid_argument("New shape must have the same total size as the current shape.");
        }

        std::copy(new_dims.begin(), new_dims.end(), dims_m.begin());

        return *this;
    }

    size_t get_total_size() const {
        return total_size();
    }

    const std::vector<T>& get_data() const {
        return data_m;
    }

    tenzor operator+= (const tenzor<T, N>& other) {
        if (this->dims_m != other.dims_m) {
            throw std::invalid_argument("Matrices must have the same dimensions");
        }

        for (size_t i = 0; i < total_size(); ++i) {
            this->data_m[i] += other.data_m[i];
        }

        return *this;
    }

    tenzor operator-= (const tenzor<T, N>& other) {
        if (this->dims_m != other.dims_m) {
            throw std::invalid_argument("Matrices must have the same dimensions");
        }

        for (size_t i = 0; i < total_size(); ++i) {
            this->data_m[i] -= other.data_m[i];
        }

        return *this;
    }

    tenzor<T, N> multiply(const tenzor<T, N>& other) const {
        auto lhs_shape = this->dims_m;
        auto rhs_shape = other.dims_m;

        if (lhs_shape[1] != rhs_shape[0]) {
            throw std::invalid_argument("Incompatible dimensions for matrix multiplication.");
        }

        size_t m = lhs_shape[0];
        size_t n = rhs_shape[1];
        size_t p = lhs_shape[1];
        std::vector<T> result_data(m * n, 0);

        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                for (size_t k = 0; k < p; ++k) {
                    result_data[i * n + j] += this->data_m[i * p + k] * other.data_m[k * n + j];
                }
            }
        }

        return tenzor<T, N>(result_data, m, n);
    }

    tenzor<T, N> invert() const {
        if (dims_m[0] != dims_m[1]) {
            throw std::invalid_argument("Matrix must be square to invert.");
        }

        size_t n = dims_m[0];
        tenzor<T, N> result{*this};
        tenzor<T, N> identity{dims_m};

        for (size_t i = 0; i < n; ++i) {
            identity(i, i) = 1;
        }

        for (size_t i = 0; i < n; ++i) {
            size_t maxRow = i;
            for (size_t k = i + 1; k < n; ++k) {
                if (std::fabs(result.data_m[k * n + i]) > std::fabs(result.data_m[maxRow * n + i])) {
                    maxRow = k;
                }
            }

            for (size_t k = 0; k < n; ++k) {
                std::swap(result.data_m[maxRow * n + k], result.data_m[i * n + k]);
                std::swap(identity.data_m[maxRow * n + k], identity.data_m[i * n + k]);
            }

            for (size_t k = i + 1; k < n; ++k) {
                T c = -result.data_m[k * n + i] / result.data_m[i * n + i];
                for (size_t j = i; j < n; ++j) {
                    result.data_m[k * n + j] += c * result.data_m[i * n + j];
                    identity.data_m[k * n + j] += c * identity.data_m[i * n + j];
                }
            }
        }

        for (int i = n - 1; i >= 0; --i) {
            for (int k = i - 1; k >= 0; --k) {
                T c = -result.data_m[k * n + i] / result.data_m[i * n + i];
                for (size_t j = 0; j < n; ++j) {
                    identity.data_m[k * n + j] += c * identity.data_m[i * n + j];
                }
            }
            T c = 1 / result.data_m[i * n + i];
            for (size_t j = 0; j < n; ++j) {
                identity.data_m[i * n + j] *= c;
            }
        }

        return identity;
    }


    std::vector<T> solve(const std::vector<T>& b) const {
        if (dims_m[0] != dims_m[1]) {
            throw std::invalid_argument("Matrix must be square to solve equations.");
        }
        if (b.size() != dims_m[0]) {
            throw std::invalid_argument("Vector size must match the number of rows in the matrix.");
        }

        tenzor<T, N> inv = this->invert();
        std::vector<T> result(b.size(), 0);

        for (size_t i = 0; i < dims_m[0]; i++) {
            for (size_t j = 0; j < dims_m[1]; j++) {
                result[i] += inv(i, j) * b[j];
            }
        }

        return result;
    }


    template <typename... Indices>
    T& operator()(Indices... indices) {
        static_assert(sizeof...(Indices) == N, "operator() requires exactly N indices.");
        return data_m[calculate_index(0, indices...)];
    }

    template <typename... Indices>
    const T& operator()(Indices... indices) const {
        return data_m[calculate_index(0, indices...)];
    }

    std::array<size_t, N> shape() const { return dims_m; }
};

template <typename T, size_t N>
tenzor<T, N> operator+ (const tenzor<T, N> &fst, const tenzor<T, N>& snd) {
    tenzor<T, N> result(fst);
    return result += snd;
}

template <typename T, size_t N>
tenzor<T, N> operator- (const tenzor<T, N> &fst, const tenzor<T, N>& snd) {
    tenzor<T, N> result(fst);
    return result -= snd;
}

#endif //LINALG_NDMATRIX_H
