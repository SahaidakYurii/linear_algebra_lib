// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#ifndef LINALG_NDMATRIX_H
#define LINALG_NDMATRIX_H

#include <iostream>
#include <vector>
#include <array>
#include <stdexcept>
#include <initializer_list>

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
    ndmatrix(std::initializer_list<size_t> dims, std::vector<T> data = {}) : data_m(data) {
    // ndmatrix(std::initializer_list<size_t> dims) : data_m(){
        if (dims.size() != N) {
            throw std::invalid_argument("Dimension size does not match template parameter N.");
        }

        std::copy(dims.begin(), dims.end(), dims_m.begin());
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
