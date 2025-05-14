# Project: linear algebra library
Authors (team): 
- Sahaidak Yurii
- Pavliuk Bohdan
- Samoilenko Marta

### Compilation

```bash
mkdir build && cd build
cmake ..
make
```

### Requirements
- tbb
- Eigen

### Usage
#### Library
Include to your project and use
#### Tests
Google tests framework, run
```bash
./tests
```
#### Benchmarks
Run one of the executables: `benchmark_runner` `benchmark_svd`, `sparseMatrix_benchmark`, 
#### Main
Used for testing, nothing there

### Structure

| class           | description                                           |
|-----------------|-------------------------------------------------------|
| Tensor<T, N>    | multidimentional matrix                               |
| Matrix<T>       | the matrix                                            |
| SparseMatrix<T> | the sparse matrix, uses unordered_map to store data   |
| SquareMatrix<T> | more functions which defined only for square matrices |
| Vector<T>       | custom vector                                         |
