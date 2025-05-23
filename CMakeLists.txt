cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME linear_algebra_lib)
project(${PROJECT_NAME} C CXX)

set(CMAKE_CXX_STANDARD 20)

option(LINALG_USE_THREADS "Enable internal linalg thread pool" OFF)

if(LINALG_USE_THREADS)
	add_compile_definitions(LINALG_USE_THREADS=1)
else()
	add_compile_definitions(LINALG_USE_THREADS=1)
endif()

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

include(FetchContent)

##########################################################
# Header-only library target
##########################################################
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE include)

if(LINALG_USE_THREADS)
	find_package(Threads REQUIRED)
	find_package(TBB REQUIRED)
	target_link_libraries(${PROJECT_NAME} INTERFACE
			Threads::Threads
			TBB::tbb)
endif()

##########################################################
# Dependencies: Eigen (header-only)
##########################################################
FetchContent_Declare(
		eigen
		GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
		GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)

#find_package(OpenCV REQUIRED COMPONENTS core imgproc highgui)

##########################################################
# Main executables
##########################################################
add_executable(main
		src/main.cpp
)
target_link_libraries(main
		PRIVATE ${PROJECT_NAME}
)
target_include_directories(main PRIVATE ${eigen_SOURCE_DIR})

##########################################################
# Benchmark executables
##########################################################


add_executable(benchmark_runner
		benchmarks/benchmark.cpp
		benchmarks/benchmark_utils.h
)
target_link_libraries(benchmark_runner
		PRIVATE ${PROJECT_NAME}
)
target_include_directories(benchmark_runner PRIVATE ${eigen_SOURCE_DIR})
target_include_directories(benchmark_runner PRIVATE benchmarks)

# EIGEN eigenvalue comparison
add_executable(eigenvalue_benchmark
		benchmarks/eigenvalue_benchmark.cpp
		benchmarks/benchmark_utils.h
)
target_link_libraries(eigenvalue_benchmark
		PRIVATE ${PROJECT_NAME}
)
target_include_directories(eigenvalue_benchmark PRIVATE ${eigen_SOURCE_DIR})
target_include_directories(benchmark_runner PRIVATE benchmarks)

# EIGEN svd comparison

add_executable(benchmark_svd
		benchmarks/bench_svd.cpp
		benchmarks/benchmark_utils.h
)

target_link_libraries(benchmark_svd
		PRIVATE
		${OpenCV_LIBS}
		Eigen3::Eigen
		${PROJECT_NAME}
)

target_include_directories(benchmark_svd PRIVATE ${OpenCV_INCLUDE_DIRS})
target_include_directories(benchmark_svd PRIVATE benchmarks)

# sparse vs matrix comparison
add_executable(sparseMatrix_benchmark
		benchmarks/sparseMatrix_benchmark.cpp
		benchmarks/benchmark_utils.h
)

target_link_libraries(sparseMatrix_benchmark
		PRIVATE ${PROJECT_NAME}
)

target_include_directories(sparseMatrix_benchmark PRIVATE benchmarks)

##########################################################
# GoogleTest: Unit test executable
##########################################################

FetchContent_Declare(
		gtest
		GIT_REPOSITORY https://github.com/google/googletest.git
		GIT_TAG release-1.12.1
)
FetchContent_MakeAvailable(gtest)

enable_testing()

add_executable(tests
		tests/test_vector.cpp
		tests/test_tenzor.cpp
		tests/test_matrix.cpp
		tests/test_squareMatrix.cpp
		tests/test_sparseMatrix.cpp
)

target_link_libraries(tests
		PRIVATE
		${PROJECT_NAME}
		gtest
		gtest_main
)

add_test(NAME unit_tests COMMAND tests)

##########################################################
# Optional: Install targets
##########################################################
install(TARGETS benchmark_runner eigenvalue_benchmark tests DESTINATION bin)
