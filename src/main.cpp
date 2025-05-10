#include <iostream>
#include "linalg.h"

int main(int argc, char* argv[]) {
    linalg::matrix<double> tall(3, 2, {1, 2, 3, 4, 5, 6});
    auto tall_pinv = tall.pinv();
    std::cout << tall.toString() << std::endl;

    std::cout << (tall * tall.pinv() * tall).toString() << std::endl;

    std::cout << (tall.pinv() * tall * tall.pinv()).toString() << std::endl;
    std::cout << tall.pinv().toString() << std::endl;

    std::cout << (tall.pinv() * tall).transpose().toString() << std::endl;
    std::cout << (tall * tall.pinv()).transpose().toString() << std::endl;

    std::cout << (tall * tall.pinv() * tall == tall) << std::endl;
    std::cout << static_cast<linalg::squareMatrix<double>>((tall * tall.pinv())).isSymmetrical() << std::endl;
}
