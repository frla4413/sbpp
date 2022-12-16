cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-fopenmp -Ofast")
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# to produce compile_commands.json:
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "${CMAKE_SOURCE_DIR}")

set(PROJECT_NAME "ins")
set(MAIN_FILE "ins/ins_main.cpp")
set(SRC_FILES "ins/src/ins.cpp" "ins/src/integrate_ins.cpp")
set(INCLUDE_DIR "ins/include")

project(${PROJECT_NAME})

add_executable(${PROJECT_NAME} ${MAIN_FILE} ${SRC_FILES})

set(SBPP_DIR "/home/fredrik/work/cpp_projects/sbpp/sbpp")
target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_DIR})
add_subdirectory(${SBPP_DIR} sbpp)

target_link_libraries(${PROJECT_NAME} PUBLIC sbpp)
