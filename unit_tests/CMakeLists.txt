cmake_minimum_required(VERSION 3.14)
project(testing)

# GoogleTest requires at least C++11
set(CMAKE_CXX_STANDARD 17)

# to produce compile_commands.json:
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "${CMAKE_SOURCE_DIR}")

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  testing
  test_main.cpp
)

set(SBPP_DIR "/home/fredrik/work/cpp_projects/sbpp/sbpp")
target_include_directories(testing PUBLIC "include/")
add_subdirectory(${SBPP_DIR} sbpp)

target_link_libraries(
  testing
  gtest_main
  sbpp
)

include(GoogleTest)
gtest_discover_tests(testing)
