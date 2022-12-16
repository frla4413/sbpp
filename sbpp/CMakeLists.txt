# set location of tqdm
set(TQDM_DIR "/home/fredrik/work/cpp_projects/external_cpp_libraries/tqdm-cpp")
set(CMAKE_CXX_STANDARD 17)

#find mkl
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

add_library(sbpp STATIC src/array.cpp
                        src/mbarray.cpp
                        src/mesh.cpp
                        src/sbp.cpp
                        src/sbp21.cpp
                        src/sbp42.cpp
                        src/mbgrid.cpp
                        src/mbsbp.cpp
                        src/basics.cpp
                        src/interp.cpp
                        src/interp21.cpp
                        src/interp42.cpp
                        src/integrate.cpp)


find_package(MKL REQUIRED)

target_include_directories(sbpp PUBLIC "-DMKL_ILP64 -m64 -I"
                                       ${MKL_INCLUDE_DIR}
                                       ${TQDM_DIR}
                                       include)

target_link_libraries(sbpp PUBLIC "-Wl,--start-group"
                                  ${MKL_LIBRARIES}
               "-Wl,--end-group -lgomp -lpthread -lm -ldl")
