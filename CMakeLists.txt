#cmake_minimum_required(VERSION 2.8.7)
cmake_minimum_required(VERSION 3.12)
project(ALGA)

#SET(CMAKE_CXX_FLAGS "-D_GLIBCXX_DEBUG -static -pthread -g -std=c++14 -O0 -Wall") # compilation for debugging
#SET(CMAKE_CXX_FLAGS " -pthread -g -std=c++14 -O0") # compilation for profiling
#SET(CMAKE_CXX_FLAGS "-static -pthread -std=c++14 -O0") # compilation for profiling
SET(CMAKE_CXX_FLAGS "-pthread -std=c++14 -lm -O3") # compilation for faster execution -O2 uses all compiler optimizations, -O0 uses no optimizations at all
#SET(CMAKE_CXX_FLAGS "-m64 -lm -O0 -g -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC -D_GLIBCXX_SANITIZE_VECTOR -D_GLIBCXX_CONCEPT_CHECKS -D_GLIBCXX_DEBUG_ASSERT -D_GLIBCXX_DEBUG_PEDASSERT -D_GLIBCXX_DEBUG_VERIFY -pedantic -march=native -ffast-math -Wextra") # compilation for faster execution -O2 uses compiler optimizations, -O0 uses no optimizations at all

#SET(CMAKE_CXX_FLAGS "-static -pthread -std=c++14 -lm -O3")

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")


include_directories(include)

file(GLOB SOURCES "src/*.cpp" "src/*/*.cpp" "src/*/*/*.cpp" "src/*/*/*/*.cpp" "src/*/*/*/*.cpp")

add_executable(ALGA ${SOURCES})

target_link_libraries(ALGA stdc++fs)