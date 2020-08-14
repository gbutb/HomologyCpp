/**
 * main_test.cpp
 * -------------
 * Contains unit tests.
 */

#include "gtest/gtest.h"

int test(int a) {
  return a;
}

TEST(testTest, Test3) {
    EXPECT_EQ(3, test(3));
}

int main(int argc, char**argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
