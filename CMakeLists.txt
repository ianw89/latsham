cmake_minimum_required(VERSION 3.16)
project(latsham CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(HDF5_ROOT "/usr/local/anaconda3n/pkgs/hdf5-1.12.1-h2b7332f_3" CACHE PATH "HDF5 root directory")
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(GSL REQUIRED)
find_package(Boost REQUIRED)
find_package(OpenMP REQUIRED)

add_library(HighFive INTERFACE)
target_include_directories(HighFive INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/extern
    ${HDF5_INCLUDE_DIRS}
)
target_link_libraries(HighFive INTERFACE ${HDF5_C_LIBRARIES})

add_executable(latsham src/main.cpp)
target_include_directories(latsham PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${Boost_INCLUDE_DIRS}
)
target_link_libraries(latsham PRIVATE HighFive GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX stdc++fs)

add_executable(tests src/tests.cpp)
target_include_directories(tests PRIVATE 
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${Boost_INCLUDE_DIRS}
)
target_link_libraries(tests PRIVATE HighFive GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX stdc++fs)
