cmake_minimum_required(VERSION 3.11)
project(darksun
	VERSION 0.1.0
	DESCRIPTION "Project for studying DM from an SU(N) gauge theory with a single quark flavor"
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(fmt REQUIRED)
find_package(Threads REQUIRED)

set(THREAD_PREFER_PTHREAD_FLAG ON)

#============================================================================
#---- Construct the library -------------------------------------------------
#============================================================================

set(DARKSUN_LIB darksun)
add_library(${DARKSUN_LIB} INTERFACE)
target_include_directories(${DARKSUN_LIB} INTERFACE  "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(${DARKSUN_LIB} PUBLIC)

set(STIFF_LIB stiff)
add_library(${STIFF_LIB} INTERFACE)
target_include_directories(${STIFF_LIB} INTERFACE  "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(${STIFF_LIB} PUBLIC)

#============================================================================
#---- Configure Google Test -------------------------------------------------
#============================================================================

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
if (result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif ()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
if (result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif ()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif ()

#============================================================================
#---- Build testing files ---------------------------------------------------
#============================================================================

set(TEST_FILES
		"test_scanner"
	"test_model"
	"test_radau")

foreach(tfile ${TEST_FILES})
	add_executable(${tfile} "test/${tfile}.cpp")
	target_link_libraries(${tfile} PUBLIC
		${DARKSUN_LIB} ${STIFF_LIB}
		gtest_main
		GSL::gsl
		GSL::gslcblas
		Eigen3::Eigen
			fmt::fmt
		${Boost_LIBRARIES})
	set_target_properties(${tfile} PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)
endforeach()

#============================================================================
#---- Build Applications ----------------------------------------------------
#============================================================================

set(APP_FILES
		"bm"
		"c"
		"c_vs_n"
		"lec1"
		"lec2"
		"mu_eta"
		"mu_del"
		"xi_inf"
	"generate_cs_data")

foreach(afile ${APP_FILES})
	add_executable(${afile} "apps/${afile}.cpp")
	target_link_libraries(${afile} PUBLIC
			${DARKSUN_LIB} 
			${STIFF_LIB}
			gtest_main
			GSL::gsl
			GSL::gslcblas
			Eigen3::Eigen
			fmt::fmt
			${Boost_LIBRARIES})
	set_target_properties(${afile} PROPERTIES
			RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)
endforeach()
