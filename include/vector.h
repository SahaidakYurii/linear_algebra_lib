// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#ifndef LINALG_VECTOR_H
#define LINALG_VECTOR_H
#include <vector>
#include <cmath>

template<typename T>
using vector = std::vector<T>;

template<typename T>
T vectorNorm(const vector<T>& v) {
    T sum = 0;
    for (auto &val : v) sum += val * val;
    return std::sqrt(sum);
}

#endif //LINALG_VECTOR_H
