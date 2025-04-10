// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H
#include <sstream>
#include "tensor.h"

template <typename T>
class matrix : public tensor<T, 2> {
protected:
    size_t rows_m, cols_m;

public:
    using tensor<T, 2>::operator();

    matrix(const size_t rows, const size_t cols, vector<T> data = {})
        : tensor<T, 2>(data, rows, cols),
          rows_m(rows), cols_m(cols) {}

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

    matrix<T> transpose() {
        matrix<T> temp(rows_m, cols_m);
        for (size_t r = 0; r < rows_m; r++) {
            for (size_t c = 0; c < cols_m; c++) {
                temp(c, r) = (*this)(r, c);
            }
        }

        return temp;
    }
};


#endif //LINALG_MATRIX_H
