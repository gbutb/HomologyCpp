cmake_minimum_required(VERSION 3.5)
project(Homology VERSION 1.0.0 LANGUAGES CXX)

# Set default options
option(BUILD_TESTS "BUILD_TESTS" ON)
option(BUILD_EXAMPLES "BUILD_EXAMPLES" ON)

IF(APPLE)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
ENDIF()

# Set to C++14
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++14 -fno-elide-constructors")

# Initialize required packages
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

# Enable debugging mode
if (DEBUG_MODE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

# Load packages
find_package(VCG REQUIRED)

# Include project
include_directories(
    ${PROJECT_SOURCE_DIR}
    ${Eigen_INCLUDE_DIRS}
    ${VCG_INCLUDE_DIRS})

if (BUILD_TESTS)
    add_subdirectory(tests/)
endif()

if (BUILD_EXAMPLES)
    add_subdirectory(examples/)
endif()
