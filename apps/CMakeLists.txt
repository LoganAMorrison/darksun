cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_STANDARD 14)

find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(EIGEN3 REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})

add_executable(scan scan.cpp)
target_link_libraries(scan
        darksun
        ${Boost_LIBRARIES}
        Eigen3::Eigen
        GSL::gsl
        GSL::gslcblas)
set_target_properties(scan PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)
