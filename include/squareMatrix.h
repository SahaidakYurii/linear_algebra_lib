//
// Created by yurii-sahaidak on 4/5/25.
//

#include "matrix.h"

#ifndef SQUAREMATRIX_H
#define SQUAREMATRIX_H

template <typename T>
class squareMatrix : public matrix<T> {
    size_t a;
protected:
    squareMatrix<T> submatrix(size_t row, size_t col) {
        std::vector<T> sub_data{};
        sub_data.reserve((this->rows_m - 1) * (this->cols_m - 1));
        for (size_t r = 0; r < this->rows_m; r++) {
            if (r == row) continue;
            for (size_t c = 0; c < this->cols_m; c++) {
                if (c == col) continue;

                sub_data.push_back(this->data_m[r * this->cols_m + c]);
            }
        }

        return squareMatrix<T>{a-1, sub_data};
    }
public:
    explicit squareMatrix(const size_t a, std::vector<T> data = {}):
        matrix<T>(a, a, data), a(a) {}

    explicit squareMatrix(const matrix<T> &other) :
        matrix<T>(other.reshape(other.cols_m < other.rows_m ? other.cols_m : other.rows_m)),
        a(other.cols_m < other.rows_m ? other.cols_m : other.rows_m)
            {}

    bool isSquare() {
        return this->rows_m == this->cols_m;
    }

    T minor(size_t row, size_t col) {
        if (!isSquare()) throw std::invalid_argument("Cofactor is only defined for square matrices");

        squareMatrix<T> sub = submatrix(row, col);
        return sub.determinant();
    }

    T cofactor(size_t row, size_t col) {
        return (((row + col) % 2 == 0) ? 1 : -1) * minor(row, col);
    }

    T determinant() {
        if (!isSquare()) throw std::invalid_argument("Determinant is only defined for square matrices");
        if (this->rows_m == 1)
            return this->data_m[0];

        if (this->rows_m == 2)
            return this->data_m[0] * this->data_m[3] - this->data_m[1] * this->data_m[2];

        T res = 0;
        for (size_t col = 0; col < this->cols_m; ++col) {
            T cof = cofactor(0, col);
            T val = this->operator()(0, col);
            res += val * cof;
        }

        return res;
    }
};

#endif //SQUAREMATRIX_H
