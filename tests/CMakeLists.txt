# Set up Google Test
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIR})
include(GoogleTest)

# Add targets
add_executable(main_test main_test.cpp)
target_link_libraries(main_test ${GTEST_LIBRARY})

# Build
gtest_add_tests(
    TARGET main_test)
