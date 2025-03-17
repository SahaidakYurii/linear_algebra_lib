#include <iostream>
#include "matrix.h"


int main(int argc, char* argv[]) {
    matrix<char> A(4, 3);
    for (size_t r = 0; r < 4; r++) {
        for (size_t c = 0; c < 3; c++) {
            A(r, c) = 'a';
        }
    }

    A.reshape(6, 2);
    std::cout << A.toString() << std::endl;
}
