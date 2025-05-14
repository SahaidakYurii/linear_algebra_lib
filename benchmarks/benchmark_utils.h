

#ifndef BENCHMARK_UTILS_H
#define BENCHMARK_UTILS_H

#include <chrono>
#include <atomic>

template<class D>
inline double to_us(const D& d)
{
    return std::chrono::duration_cast<std::chrono::microseconds>(d).count();
}

inline std::chrono::high_resolution_clock::time_point get_current_time_fenced()
{
    std::atomic_thread_fence(std::memory_order_seq_cst);
    auto res_time = std::chrono::high_resolution_clock::now();
    std::atomic_thread_fence(std::memory_order_seq_cst);
    return res_time;
}

template<typename Func>
double benchmark(Func func, int repetitions = 10) {
    double total_time = 0.0;
    for (int i = 0; i < repetitions; ++i) {
        auto start = get_current_time_fenced();
        func();
        auto end = get_current_time_fenced();
        total_time += to_us(end - start);
    }
    return total_time / repetitions;
}

#endif //BENCHMARK_UTILS_H
