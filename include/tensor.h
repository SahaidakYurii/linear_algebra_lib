// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#ifndef LINALG_NDMATRIX_H
#define LINALG_NDMATRIX_H

#include <iostream>
#include <array>
#include <stdexcept>
#include <initializer_list>
#include <cmath>

#include "vector.h"

template <typename T, size_t N>
class tensor {
protected:
    vector<T> data_m;
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
            throw std::invalid_argument("Number of arguments does not match dimentionality of the tensor");
        }
        return 0;
    }

public:
    template <typename... Indices>
    explicit tensor(Indices ... indices) : data_m() {
        static_assert(sizeof...(Indices) == N, "Number of indices must match the dimension N");
        static_assert((std::is_convertible_v<Indices, size_t> && ...), "All indices must be convertible to size_t");

        dims_m = { static_cast<size_t>(indices)... };
        data_m.resize(total_size());
    }

    template <typename... Indices>
    explicit tensor(const vector<T>& data, Indices... indices)
            : tensor(indices...){
        data_m = data;
        data_m.resize(total_size());

    }

    tensor<T, N>& reshape(std::initializer_list<size_t> new_dims) {
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

    const vector<T>& get_data() const {
        return data_m;
    }

    bool operator==(const tensor<T, N>& other) const {
        for (size_t i = 0; i < data_m.size(); ++i) {
            if (data_m[i] != other.data_m[i]) {
                return false;
            }
        }
        return true;
    }

    vector<T>& get_data() {
        return data_m;
    }

    tensor operator+= (const tensor<T, N>& other) {
        if (this->dims_m != other.dims_m) {
            throw std::invalid_argument("Matrices must have the same dimensions");
        }

        for (size_t i = 0; i < total_size(); ++i) {
            this->data_m[i] += other.data_m[i];
        }

        return *this;
    }

    tensor operator-= (const tensor<T, N>& other) {
        if (this->dims_m != other.dims_m) {
            throw std::invalid_argument("Matrices must have the same dimensions");
        }

        for (size_t i = 0; i < total_size(); ++i) {
            this->data_m[i] -= other.data_m[i];
        }

        return *this;
    }

    tensor<T, N> multiply(const tensor<T, N>& other) const {
        auto lhs_shape = this->dims_m;
        auto rhs_shape = other.dims_m;

        if (lhs_shape[1] != rhs_shape[0]) {
            throw std::invalid_argument("Incompatible dimensions for matrix multiplication.");
        }

        size_t m = lhs_shape[0];
        size_t n = rhs_shape[1];
        size_t p = lhs_shape[1];
        vector<T> result_data(m * n, 0);

        for (size_t i = 0; i < m; ++i) {
            for (size_t j = 0; j < n; ++j) {
                for (size_t k = 0; k < p; ++k) {
                    result_data[i * n + j] += this->data_m[i * p + k] * other.data_m[k * n + j];
                }
            }
        }

        return tensor<T, N>(result_data, m, n);
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
tensor<T, N> operator+ (const tensor<T, N> &fst, const tensor<T, N>& snd) {
    tensor<T, N> result(fst);
    return result += snd;
}

template <typename T, size_t N>
tensor<T, N> operator- (const tensor<T, N> &fst, const tensor<T, N>& snd) {
    tensor<T, N> result(fst);
    return result -= snd;
}

#endif //LINALG_NDMATRIX_H
