// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H

#include "ndmatrix.h"
#include <vector>
#include <sstream>

template <typename T>
class matrix : ndmatrix<T, 2> {
    size_t rows_m{}, cols_m{};
    std::vector<T> data_m;

public:
    matrix(const size_t rows, const size_t cols, std::vector<T> data = {})
        : rows_m(rows), cols_m(cols), data_m(data) {
        data_m.resize(rows_m * cols_m);
    };
    matrix() = default;
    ~matrix() = default;

    T& operator()(const size_t row, const size_t col) {
        return data_m[row * cols_m + col];
    }
    const T& operator()(const size_t row, const size_t col) const {
        return data_m[row * cols_m + col];
    }

    void reshape(const size_t rows, const size_t cols) {
        if (rows == rows_m || cols == cols_m) {
            return;
        }

        auto temp_data = std::vector<T>(rows * cols);

        // copies the common part, if 3x4 reshaped to 2x6, the 2x4 part will remain same
        for (size_t r = 0; r < (rows < rows_m ? rows : rows_m); r++) {
            for (size_t c = 0; c < (cols < cols_m ? cols : cols_m); c++) {
                temp_data[r * cols + c] = data_m[r * cols_m + c];
            }
        }
        rows_m = rows; cols_m = cols;
        data_m = temp_data;
    }

    std::string toString() {
        std::stringstream ss;
        for (size_t r = 0; r < rows_m; r++) {
            for (size_t c = 0; c < cols_m; c++) {
                ss << data_m[r * cols_m + c] << '\t';
            }
            ss << '\n';
        }

        return ss.str();
    }
};


#endif //LINALG_MATRIX_H
