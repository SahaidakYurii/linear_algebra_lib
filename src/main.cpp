#include <iostream>
#include "linalg.h"

int main(int argc, char* argv[]) {
    linalg::vector<linalg::vector<int>> jagged = {
        {1, 2},
        {3},
        {4, 5, 6}
    };

    linalg::matrix<int> m(jagged);

    std::cout << m.toString() << std::endl;
}
