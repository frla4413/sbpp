cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

#set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# to produce compile_commands.json:
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "${CMAKE_SOURCE_DIR}")

set(PROJECT_NAME "test")
project(${PROJECT_NAME} LANGUAGES CXX C)

set(MAIN_FILE "main.cpp")

project(${PROJECT_NAME})

add_executable(${PROJECT_NAME} ${MAIN_FILE})

set(SBPP_DIR "/home/fredrik/work/cpp_projects/sbpp/sbpp")
target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_DIR})
add_subdirectory(${SBPP_DIR} sbpp)

target_link_libraries(${PROJECT_NAME} PUBLIC sbpp)
