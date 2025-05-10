// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#include <gtest/gtest.h>
#include "linalg.h"
#include <stdexcept>
#include <initializer_list>
#include <iterator>
#include <algorithm>
#include <string>

TEST(vector, DefaultConstructor) {
    linalg::vector<int> v;
    EXPECT_TRUE(v.is_empty());
    EXPECT_EQ(v.size(), 0);
}

TEST(vector, FillConstructor) {
    linalg::vector<int> v(static_cast<size_t>(5), 42);
    EXPECT_EQ(v.size(), 5);
    for (int i = 0; i < 5; ++i)
        EXPECT_EQ(v[i], 42);
}

TEST(vector, IteratorConstructor) {
    std::vector<int> source = {1, 2, 3, 4};
    linalg::vector<int> v(source.begin(), source.end());
    EXPECT_EQ(v.size(), 4);
    EXPECT_EQ(v[2], 3);
}

TEST(vector, InitializerListConstructor) {
    linalg::vector<std::string> v{"a", "b", "c"};
    EXPECT_EQ(v.size(), 3);
    EXPECT_EQ(v[0], "a");
    EXPECT_EQ(v[1], "b");
    EXPECT_EQ(v[2], "c");
}

TEST(vector, CopyConstructor) {
    linalg::vector<int> v1{1, 2, 3};
    linalg::vector<int> v2 = v1;
    EXPECT_EQ(v2, v1);
}

TEST(vector, MoveConstructor) {
    linalg::vector<int> temp{5, 6};
    linalg::vector<int> v = std::move(temp);
    EXPECT_EQ(v.size(), 2);
    EXPECT_EQ(v[0], 5);
    EXPECT_EQ(v[1], 6);
}

TEST(vector, CopyAssignment) {
    linalg::vector<int> a{1, 2, 3};
    linalg::vector<int> b;
    b = a;
    EXPECT_EQ(a, b);
}

TEST(vector, MoveAssignment) {
    linalg::vector<int> a{1, 2};
    linalg::vector<int> b;
    b = std::move(a);
    EXPECT_EQ(b.size(), 2);
    EXPECT_EQ(b[0], 1);
    EXPECT_EQ(b[1], 2);
}

TEST(vector, IndexOperator) {
    linalg::vector<char> v{'a', 'b'};
    EXPECT_EQ(v[1], 'b');
}

TEST(vector, AtThrowsOutOfBounds) {
    linalg::vector<int> v{1, 2, 3};
    EXPECT_THROW(v.at(10), std::out_of_range);
}

TEST(vector, FrontBack) {
    linalg::vector<int> v{10, 20, 30};
    EXPECT_EQ(v.front(), 10);
    EXPECT_EQ(v.back(), 30);
}

TEST(vector, Iterators) {
    linalg::vector<int> v{1, 2, 3};
    auto it = std::find(v.begin(), v.end(), 2);
    EXPECT_NE(it, v.end());
    EXPECT_EQ(*it, 2);
}

TEST(vector, BackInserterCompatibility) {
    linalg::vector<int> v;
    std::vector<int> source = {1, 2, 3};
    std::copy(source.begin(), source.end(), std::back_inserter(v));
    EXPECT_EQ(v.size(), 3);
    EXPECT_EQ(v.back(), 3);
}

TEST(vector, CapacityManagement) {
    linalg::vector<int> v;
    v.reserve(10);
    EXPECT_GE(v.capacity(), 10);
    v.shrink_to_fit();
    EXPECT_LE(v.capacity(), 10);
}

TEST(vector, ResizeAndClear) {
    linalg::vector<int> v{1, 2, 3};
    v.resize(5);
    EXPECT_EQ(v.size(), 5);
    v.clear();
    EXPECT_TRUE(v.is_empty());
}

TEST(vector, ResizeShrink) {
    linalg::vector<int> v{1, 2, 3, 4, 5};
    v.resize(3);
    EXPECT_EQ(v.size(), 3);
    EXPECT_EQ(v[0], 1);
    EXPECT_EQ(v[2], 3);
}

TEST(vector, ResizeGrowDefault) {
    linalg::vector<int> v{1, 2};
    v.resize(5);
    EXPECT_EQ(v.size(), 5);
    EXPECT_EQ(v[2], {});
}

TEST(vector, InsertSingle) {
    linalg::vector<int> v{1, 3};
    auto it = v.insert(v.begin() + 1, 2);
    EXPECT_EQ(*it, 2);
    EXPECT_EQ(v[1], 2);
}

TEST(vector, InsertAtBeginEmpty) {
    linalg::vector<int> v;
    v.insert(v.begin(), 42);
    EXPECT_EQ(v.size(), 1);
    EXPECT_EQ(v[0], 42);
}

TEST(vector, InsertAtEnd) {
    linalg::vector<int> v{1, 2, 3};
    v.insert(v.end(), 4);
    EXPECT_EQ(v.size(), 4);
    EXPECT_EQ(v[3], 4);
}

TEST(vector, InsertRange) {
    linalg::vector<int> v{1, 5};
    std::vector<int> to_insert{2, 3, 4};
    v.insert(v.begin() + 1, to_insert.begin(), to_insert.end());
    EXPECT_EQ(v[2], 3);
}

TEST(vector, EraseSingle) {
    linalg::vector<int> v{1, 2, 3};
    auto it = v.erase(v.begin() + 1);
    EXPECT_EQ(*it, 3);
    EXPECT_EQ(v.size(), 2);
}

TEST(vector, EraseRange) {
    linalg::vector<int> v{1, 2, 3, 4, 5};
    auto it = v.erase(v.begin() + 1, v.begin() + 4);
    EXPECT_EQ(*it, 5);
    EXPECT_EQ(v.size(), 2);
}

TEST(vector, EraseFromEmpty) {
    linalg::vector<int> v;
    EXPECT_THROW(v.erase(v.begin()), std::out_of_range);
}

TEST(vector, PushPopEmplaceBack) {
    linalg::vector<std::string> v;
    v.push_back("hi");
    v.emplace_back("there");
    EXPECT_EQ(v.size(), 2);
    EXPECT_EQ(v.back(), "there");
    v.pop_back();
    EXPECT_EQ(v.size(), 1);
}

TEST(vector, ComparisonOperators) {
    linalg::vector<int> a{1, 2, 3};
    linalg::vector<int> b{1, 2, 4};
    EXPECT_TRUE(a < b);
    EXPECT_TRUE(a != b);
    EXPECT_FALSE(a == b);
}

TEST(vector, SwapMethod) {
    linalg::vector<int> a{1, 2}, b{3, 4};
    a.swap(b);
    EXPECT_EQ(a[0], 3);
    EXPECT_EQ(b[0], 1);
}

TEST(vector, SwapWithEmptyVector) {
    linalg::vector<int> v1{1, 2, 3};
    linalg::vector<int> v2;
    v1.swap(v2);
    EXPECT_EQ(v1.size(), 0);
    EXPECT_EQ(v2.size(), 3);
    EXPECT_EQ(v2[0], 1);
}

TEST(vector, NestedVectorsBasic) {
    linalg::vector<linalg::vector<int>> vv;
    linalg::vector<int> inner1 = {1, 2, 3};
    linalg::vector<int> inner2 = {4, 5};

    vv.push_back(inner1);
    vv.push_back(inner2);

    ASSERT_EQ(vv.size(), 2);
    EXPECT_EQ(vv[0].size(), 3);
    EXPECT_EQ(vv[1][1], 5);
}

TEST(vector, NestedVectorsEmplaceBack) {
    linalg::vector<linalg::vector<int>> vv;
    vv.emplace_back(linalg::vector<int>{1, 2});
    vv.emplace_back(linalg::vector<int>{3, 4});

    EXPECT_EQ(vv.size(), 2);
    EXPECT_EQ(vv[0].front(), 1);
    EXPECT_EQ(vv[1].back(), 4);
}

TEST(vector, PushBackSelfBack) {
    linalg::vector<int> v;
    for (int i = 0; i < 10; ++i) {
        v.push_back(i);
    }

    int last = v.back();
    v.push_back(last);

    EXPECT_EQ(v.size(), 11);
    EXPECT_EQ(v.back(), 9);
}

TEST(vector, NestedPushBackSelfBack) {
    linalg::vector<linalg::vector<int>> v;
    linalg::vector<int> tmp = {1, 2};
    v.push_back(tmp);
    v.push_back(v.back()); // Copy of the last vector

    ASSERT_EQ(v.size(), 2);
    EXPECT_EQ(v[1], v[0]); // Must be deep copy
}

