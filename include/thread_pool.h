#ifndef LINALG_THREAD_POOL_H
#define LINALG_THREAD_POOL_H

#include <vector>
#include <queue>
#include <thread>
#include <future>
#include <functional>
#include <condition_variable>

#if !defined(LINALG_USE_THREADS)
#  define LINALG_USE_THREADS 0   // default = ON
#endif

namespace linalg::detail {

    class ThreadPool {
    public:
        static ThreadPool& instance()
        {
            static ThreadPool pool;
            return pool;
        }

        template<class F>
        auto enqueue(F&& f) -> std::future<decltype(f())>
        {
#if LINALG_USE_THREADS
            using R = decltype(f());
            auto task = std::make_shared<std::packaged_task<R()>>(std::forward<F>(f));
            std::future<R> fut = task->get_future();

            {
                std::unique_lock lock(mutex_);
                jobs_.emplace([task](){ (*task)(); });
            }
            cv_.notify_one();
            return fut;
#else
            std::promise<decltype(f())> p;
            p.set_value(f());
            return p.get_future();
#endif
        }

        size_t size() const noexcept { return threads_.size(); }

    private:
        ThreadPool()
        {
    #if LINALG_USE_THREADS
            const unsigned n = std::max(1u, std::thread::hardware_concurrency());
            for (unsigned i = 0; i < n; ++i)
                threads_.emplace_back([this] { worker(); });
    #endif
        }
        ~ThreadPool()
        {
    #if LINALG_USE_THREADS
            {   std::unique_lock lock(mutex_);
                done_ = true;
            }
            cv_.notify_all();
            for (auto& t : threads_) t.join();
    #endif
        }

        void worker()
        {
            while (true)
            {
                std::function<void()> job;
                {   std::unique_lock lock(mutex_);
                    cv_.wait(lock, [this]{ return done_ || !jobs_.empty(); });
                    if (done_ && jobs_.empty()) return;
                    job = std::move(jobs_.front());
                    jobs_.pop();
                }
                job();
            }
        }

        bool done_ = false;
        std::vector<std::thread>            threads_;
        std::queue<std::function<void()>>   jobs_;
        std::mutex                          mutex_;
        std::condition_variable             cv_;
    };

}
#endif // LINALG_THREAD_POOL_H
