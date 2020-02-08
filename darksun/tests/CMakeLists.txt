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

find_package(Boost REQUIRED)
find_package(GSL REQUIRED)
find_package(EIGEN3 REQUIRED)

include_directories(
        ${Boost_INCLUDE_DIRS}
)

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(test_thermodynamic_particle test_thermodynamic_particle.cpp)
target_link_libraries(test_thermodynamic_particle
        darksun
        gtest_main
        Eigen3::Eigen
        ${Boost_LIBRARIES}
        GSL::gsl
        GSL::gslcblas)

set_target_properties(test_thermodynamic_particle PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)


# Now simply link against gtest or gtest_main as needed. Eg
add_executable(test_standard_model test_standard_model.cpp)
target_link_libraries(test_standard_model
        gtest_main
        darksun
        Eigen3::Eigen
        ${Boost_LIBRARIES}
        GSL::gsl
        GSL::gslcblas)

set_target_properties(test_standard_model PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)


# Now simply link against gtest or gtest_main as needed. Eg
add_executable(test_compute_xi test_compute_xi.cpp)
target_link_libraries(test_compute_xi
        darksun
        gtest_main
        ${Boost_LIBRARIES})

set_target_properties(test_compute_xi PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)

# Now simply link against gtest or gtest_main as needed. Eg
add_executable(test_thermal_cross_sections test_thermal_cross_sections.cpp)
target_link_libraries(test_thermal_cross_sections
        gtest_main
        darksun
        Eigen3::Eigen
        ${Boost_LIBRARIES}
        GSL::gsl
        GSL::gslcblas)

set_target_properties(test_thermal_cross_sections PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)


# Now simply link against gtest or gtest_main as needed. Eg
add_executable(test_odeint test_odeint.cpp)
target_link_libraries(test_odeint
        gtest_main
        darksun
        Eigen3::Eigen
        ${Boost_LIBRARIES}
        GSL::gsl
        GSL::gslcblas)


add_executable(test_models test_models.cpp)
target_link_libraries(test_models
        gtest_main
        darksun
        Eigen3::Eigen
        ${Boost_LIBRARIES}
        GSL::gsl
        GSL::gslcblas)

set_target_properties(test_models PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}../../bin)