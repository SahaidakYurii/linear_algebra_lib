cmake_minimum_required(VERSION 3.15)

#! CHANGE YOUR PROJECT NAME
set(PROJECT_NAME linear_algebra_lib)
project(${PROJECT_NAME} C CXX)

set(CMAKE_CXX_STANDARD 20)

# Default to Release build type
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

##########################################################
# Project files, libraries, and configurations
##########################################################

# Add Eigen (header-only)
include(FetchContent)
FetchContent_Declare(
		eigen
		GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
		GIT_TAG 3.4.0
)
FetchContent_MakeAvailable(eigen)

# --------------------------------------------
# Main executable
# --------------------------------------------
add_executable(${PROJECT_NAME}
		src/main.cpp
		include/linalg.h
		include/matrix.h
		include/tensor.h
		include/squareMatrix.h
		include/vector.h
)

target_include_directories(${PROJECT_NAME} PRIVATE include)

# --------------------------------------------
# Benchmark executable
# --------------------------------------------
add_executable(benchmark_runner
		benchmarks/benchmark.cpp
		include/linalg.h
		include/matrix.h
		include/tensor.h
		include/vector.h
)

target_include_directories(benchmark_runner
		PRIVATE include
		PRIVATE ${eigen_SOURCE_DIR}
)

add_executable(eigenvalue_benchmark
		benchmarks/eigenvalue_benchmark.cpp
		include/squareMatrix.h
		include/matrix.h
		include/tensor.h
		include/vector.h
)

target_include_directories(eigenvalue_benchmark PRIVATE include ${eigen_SOURCE_DIR})

# --------------------------------------------
# Test executable (NO GoogleTest used)
# --------------------------------------------
add_executable(tests
		tests/test_main.cpp
)

target_include_directories(tests PRIVATE include)

# --------------------------------------------
# Install targets
# --------------------------------------------
INSTALL(PROGRAMS
		$<TARGET_FILE:${PROJECT_NAME}>
		$<TARGET_FILE:benchmark_runner>
		DESTINATION bin)
