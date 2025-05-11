// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: https://pvs-studio.com

#ifndef LINALG_VECTOR_H
#define LINALG_VECTOR_H

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iterator>
#include <initializer_list>
#include <compare>

namespace linalg {
    template<typename T>
    class vector {
    private:
        T* data_;
        size_t size_;
        size_t capacity_;

    public:
        using value_type = T;

        vector() noexcept
            : data_(nullptr), size_(0), capacity_(0) {}

        vector(size_t count, const T& value)
            : data_(nullptr), size_(0), capacity_(0) {
            reserve(count);
            for (size_t i = 0; i < count; ++i)
                new (&data_[i]) T(value);
            size_ = count;
        }

        explicit vector(size_t len)
        : data_(nullptr), size_(0), capacity_(0) {
            reserve(len);
            for (size_t i = 0; i < len; ++i)
                new (&data_[i]) T();
            size_ = len;
        }

        template<typename InputIt>
        vector(InputIt first, InputIt last)
            : data_(nullptr), size_(0), capacity_(0) {
            size_t count = std::distance(first, last);
            reserve(count);
            std::copy(first, last, data_);
            size_ = count;
        }

        vector(std::initializer_list<T> init)
        : data_(static_cast<T*>(::operator new(init.size() * sizeof(T)))),
          size_(init.size()),
          capacity_(init.size())
        {
            size_t i = 0;
            try {
                for (const auto& item : init) {
                    new (&data_[i++]) T(item);
                }
            } catch (...) {
                for (size_t j = 0; j < i; ++j) {
                    data_[j].~T();
                }
                ::operator delete(data_);
                throw;
            }
        }

        vector(const vector& other)
        : data_(static_cast<T*>(::operator new(other.capacity_ * sizeof(T)))),
          size_(other.size_),
          capacity_(other.capacity_)
        {
            for (size_t i = 0; i < size_; ++i) {
                new (&data_[i]) T(other.data_[i]);  // deep copy each element
            }
        }

        vector(vector&& other) noexcept
            : data_(other.data_),
              size_(other.size_),
              capacity_(other.capacity_)
        {
            other.data_ = nullptr;
            other.size_ = 0;
            other.capacity_ = 0;
        }

        ~vector() {
            clear();  // Explicitly destroys all elements
            ::operator delete(data_);  // Raw memory deallocation without calling destructors
        }

        vector& operator=(const vector& other) {
            if (this != &other) {
                vector temp(other);
                swap(temp);
            }
            return *this;
        }

        vector& operator=(vector&& other) noexcept {
            if (this != &other) {
                clear();
                ::operator delete(data_);
                data_ = other.data_;
                size_ = other.size_;
                capacity_ = other.capacity_;
                other.data_ = nullptr;
                other.size_ = other.capacity_ = 0;
            }
            return *this;
        }

        T& operator[](size_t index) noexcept { return data_[index]; }
        const T& operator[](size_t index) const noexcept { return data_[index]; }

        T& at(size_t index) {
            if (index >= size_) throw std::out_of_range("vector::at");
            return data_[index];
        }
        const T& at(size_t index) const {
            if (index >= size_) throw std::out_of_range("vector::at");
            return data_[index];
        }

        T& front() {
            if (is_empty()) throw std::out_of_range("vector::front");
            return data_[0];
        }
        const T& front() const {
            if (is_empty()) throw std::out_of_range("vector::front");
            return data_[0];
        }

        T& back() {
            if (is_empty()) throw std::out_of_range("vector::back");
            return data_[size_ - 1];
        }
        const T& back() const {
            if (is_empty()) throw std::out_of_range("vector::back");
            return data_[size_ - 1];
        }

        T* begin() noexcept { return data_; }
        T* end() noexcept { return data_ + size_; }
        const T* begin() const noexcept { return data_; }
        const T* end() const noexcept { return data_ + size_; }
        const T* cbegin() const noexcept { return data_; }
        const T* cend() const noexcept { return data_ + size_; }

        auto rbegin() noexcept { return std::reverse_iterator<T*>(end()); }
        auto rend() noexcept { return std::reverse_iterator<T*>(begin()); }
        auto rbegin() const noexcept { return std::reverse_iterator<const T*>(end()); }
        auto rend() const noexcept { return std::reverse_iterator<const T*>(begin()); }
        auto rcbegin() const noexcept { return rbegin(); }
        auto rcend() const noexcept { return rend(); }

        bool is_empty() const noexcept { return size_ == 0; }
        size_t size() const noexcept { return size_; }
        size_t capacity() const noexcept { return capacity_; }

        void reserve(size_t new_cap) {
            if (new_cap <= capacity_) return;

            T* new_data = static_cast<T*>(::operator new(new_cap * sizeof(T)));

            for (size_t i = 0; i < size_; ++i) {
                new (&new_data[i]) T(std::move_if_noexcept(data_[i]));
                data_[i].~T();
            }

            ::operator delete(data_);
            data_ = new_data;
            capacity_ = new_cap;
        }

        void shrink_to_fit() {
            if (size_ < capacity_) {
                T* new_data = nullptr;
                if (size_) {
                    new_data = static_cast<T*>(::operator new(size_ * sizeof(T)));
                    for (size_t i = 0; i < size_; ++i) {
                        new (&new_data[i]) T(std::move_if_noexcept(data_[i]));
                        data_[i].~T();
                    }
                }
                ::operator delete(data_);
                data_ = new_data;
                capacity_ = size_;
            }
        }

        void clear() noexcept {
            for (size_t i = 0; i < size_; ++i)
                data_[i].~T();
            size_ = 0;
        }

        void resize(size_t count, const T& value = T()) {
            if (count < size_) {
                for (size_t i = count; i < size_; ++i)
                    data_[i].~T();
            } else if (count > size_) {
                reserve(count);
                for (size_t i = size_; i < count; ++i)
                    new (&data_[i]) T(value);
            }
            size_ = count;
        }

        void swap(vector& other) noexcept {
            std::swap(data_, other.data_);
            std::swap(size_, other.size_);
            std::swap(capacity_, other.capacity_);
        }

        T* insert(T* pos, const T& value) {
            size_t idx = pos - data_;
            if (pos < data_ || pos > end()) throw std::out_of_range("insert");
            if (size_ == capacity_) reserve(capacity_ ? capacity_ * 2 : 1);
            pos = data_ + idx;
            for (size_t i = size_; i > idx; --i) {
                new (&data_[i]) T(std::move_if_noexcept(data_[i - 1]));
                data_[i - 1].~T();
            }
            new (&data_[idx]) T(value);
            ++size_;
            return data_ + idx;
        }

        template<typename InputIt>
        T* insert(T* pos, InputIt first, InputIt last) {
            size_t idx = pos - data_;
            size_t count = std::distance(first, last);
            if (pos < data_ || pos > end()) throw std::out_of_range("insert range");
            if (size_ + count > capacity_) reserve(std::max(capacity_ * 2, size_ + count));
            pos = data_ + idx;
            for (size_t i = size_; i > idx; --i) {
                new (&data_[i + count - 1]) T(std::move_if_noexcept(data_[i - 1]));
                data_[i - 1].~T();
            }
            for (size_t i = 0; i < count; ++i)
                new (&data_[idx + i]) T(*std::next(first, i));
            size_ += count;
            return data_ + idx;
        }

        T* erase(T* pos) {
            if (pos < data_ || pos >= end()) throw std::out_of_range("erase");
            size_t idx = pos - data_;
            data_[idx].~T();
            for (size_t i = idx; i < size_ - 1; ++i) {
                new (&data_[i]) T(std::move_if_noexcept(data_[i + 1]));
                data_[i + 1].~T();
            }
            --size_;
            return data_ + idx;
        }

        T* erase(T* first, T* last) {
            if (first < data_ || last > end() || first > last) throw std::out_of_range("erase range");
            size_t idx = first - data_;
            size_t count = last - first;
            for (size_t i = 0; i < count; ++i)
                data_[idx + i].~T();
            for (size_t i = idx + count; i < size_; ++i) {
                new (&data_[i - count]) T(std::move_if_noexcept(data_[i]));
                data_[i].~T();
            }
            size_ -= count;
            return data_ + idx;
        }

        void push_back(const T& value) {
            if (size_ == capacity_) {
                // handle self reference
                T value_copy = value;
                reserve(capacity_ ? capacity_ * 2 : 1);
                new (&data_[size_++]) T(std::move(value_copy));
            } else {
                new (&data_[size_++]) T(value);
            }
        }

        void concat(const vector<T>& v2) {
            reserve(this->size() + v2.size());
            std::copy(v2.begin(), v2.end(), std::back_inserter(*this));
        }

        void pop_back() {
            if (is_empty()) throw std::out_of_range("pop_back");
            data_[--size_ - 1].~T();
        }

        template<typename... Args>
        void emplace_back(Args&&... args) {
            if (size_ == capacity_)
                reserve(capacity_ ? capacity_ * 2 : 1);
            new (&data_[size_++]) T(std::forward<Args>(args)...);
        }

        auto operator<=>(const vector& other) const {
            return std::lexicographical_compare_three_way(
                begin(), end(), other.begin(), other.end());
        }

        bool operator==(const vector& other) const {
            if (size_ != other.size_) return false;
            for (size_t i = 0; i < size_; ++i)
                if (!(data_[i] == other.data_[i]))
                    return false;
            return true;
        }
    };

    template<typename T>
    T vectorNorm(const vector<T>& v) {
        T sum = 0;
        for (auto &val : v) sum += val * val;
        return std::sqrt(sum);
    }
} // namespace linalg

#endif //LINALG_VECTOR_H
