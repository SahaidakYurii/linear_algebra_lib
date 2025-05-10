#include <iostream>
#include "linalg.h"

int main(int argc, char* argv[]) {
    // matrix<char> A(4, 3);
    // for (size_t r = 0; r < 4; r++) {
    //     for (size_t c = 0; c < 3; c++) {
    //         A(r, c) = 'a';
    //     }
    // }
    //
    // A.reshape(6, 2);
    // std::cout << A.toString() << std::endl;
    //
    // tensor<char, 2> B{{2, 2}, {'a', 'b'}};
    // for (size_t i = 0; i < B.shape()[0]; i++) {
    //     for (size_t j = 0; j < B.shape()[1]; j++) {
    //         std::cout << B(i, j) << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << std::endl;


    // tensor<int, 3> С{4, 3, 3};
    //
    // size_t i = 0;
    // for (size_t r = 0; r < 4; r++) {
    //     for (size_t c = 0; c < 3; c++) {
    //         for (size_t k = 0; k < 3; k++) {
    //             С(r, c, k) = i++;
    //         }
    //     }
    // }
    //
    // for (size_t r = 0; r < 4; r++) {
    //     for (size_t c = 0; c < 3; c++) {
    //         for (size_t k = 0; k < 3; k++) {
    //             std::cout << С(r, c, k) << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << "---------------------" << std::endl;
    //     std::cout << std::endl;
    // }
    //  Type---    --Num of dim        -------- size of each dimension
    //        |   |          ---------|
    //        V   V          V  V  V
    linalg::tensor<int, 3> m{2, 3, 4};

    size_t i = 0;
    for (size_t r = 0; r < 2; r++) {
        for (size_t c = 0; c < 3; c++) {
            for (size_t k = 0; k < 4; k++) {
                m(r, c, k) = i++;
            }
        }
    }

    for (size_t r = 0; r < 2; r++) {
        for (size_t c = 0; c < 3; c++) {
            for (size_t k = 0; k < 4; k++) {
                std::cout << m(r, c, k) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "---------------------" << std::endl;
        std::cout << std::endl;
    }

    std::cout << "-----------------------------------------" << std::endl;
    m.reshape({4, 3, 2});

    for (size_t r = 0; r < 4; r++) {
        for (size_t c = 0; c < 3; c++) {
            for (size_t k = 0; k < 2; k++) {
                std::cout << m(r, c, k) << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "---------------------" << std::endl;
        std::cout << std::endl;
    }

    linalg::squareMatrix<int> square{3, {1, 0, 0, 0, 1, 0, 0, 0, 1}};

    square.transpose();

    std::cout << square.determinant() << std::endl;
}
