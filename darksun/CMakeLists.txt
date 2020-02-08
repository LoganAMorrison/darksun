cmake_minimum_required(VERSION 3.11...3.16)

project(DarkSUN
        VERSION 0.1
        DESCRIPTION "Dark SU(N) DM theory"
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

project(DarkSUN)

find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenMP)
find_package(Eigen3 REQUIRED)

set(CMAKE_MACOSX_RPATH ON)

file(GLOB DARKSUN_SOURCES "${CMAKE_SOURCE_DIR}/darksun/src/*.cpp")

add_library(darksun ${DARKSUN_SOURCES})

target_include_directories(darksun PUBLIC "${CMAKE_SOURCE_DIR}/darksun/include")

target_link_libraries(darksun
        ${Boost_LIBRARIES}
        Eigen3::Eigen
        GSL::gsl
        GSL::gslcblas)


add_subdirectory(tests)