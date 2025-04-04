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
class ndmatrix {
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
        size_t multiplier = 1;
        for (size_t i = 0; i < iter; ++i) {
            multiplier *= dims_m[i];
        }
        iter++;
        return idx * multiplier + calculate_index(iter, indices...);
    }
    size_t calculate_index(size_t iter) const {
        if (iter != N) {
            throw std::invalid_argument("Number of arguments does not match dimentionality of the ndmatrix");
        }
        return 0;
    }



public:



    ndmatrix(const std::array<size_t, N>& dims, const std::vector<T>& data = {})
            : dims_m(dims), data_m(data) {
        data_m.resize(total_size());
    }

    // Constructor that allows initialization from an initializer list
    ndmatrix(std::initializer_list<size_t> dims_list, std::vector<T> data = {})
            : data_m(data) {
        if (dims_list.size() != N) {
            throw std::invalid_argument("Initializer list size must match template parameter N.");
        }
        std::copy(dims_list.begin(), dims_list.end(), dims_m.begin());
        data_m.resize(total_size());
    }

    ndmatrix<T, N>& reshape(std::initializer_list<size_t> new_dims) {
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

    ndmatrix<T, N> add(const ndmatrix<T, N>& other) const {
        if (this->dims_m != other.dims_m) {
            throw std::invalid_argument("Matrices must have the same dimensions for addition.");
        }

        std::vector<T> result_data(this->data_m);
        for (size_t i = 0; i < total_size(); ++i) {
            result_data[i] += other.data_m[i];
        }

        return ndmatrix<T, N>(this->dims_m, result_data);
    }

    ndmatrix<T, N> multiply(const ndmatrix<T, N>& other) const {
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

        return ndmatrix<T, N>({m, n}, result_data);
    }

    ndmatrix<T, N> invert() const {
        if (dims_m[0] != dims_m[1]) {
            throw std::invalid_argument("Matrix must be square to invert.");
        }

        size_t n = dims_m[0];
        ndmatrix<T, N> result(*this);
        ndmatrix<T, N> identity(dims_m);

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

        ndmatrix<T, N> inv = this->invert();
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
        return data_m[calculate_index(0, indices...)];
    }

    template <typename... Indices>
    const T& operator()(Indices... indices) const {
        return data_m[calculate_index(0, indices...)];
    }



    std::array<size_t, N> shape() const { return dims_m; }
};




#endif //LINALG_NDMATRIX_H
