cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
cmake_policy(VERSION 3.0.2)

project(epsilon VERSION 0.1 LANGUAGES CXX)

# Control whether to add unit test target or not
option(EPSILON_ENABLE_UNIT_TESTS "Create unit test target" ON)

# Basic library target (header only)
add_library(epsilon INTERFACE)
target_include_directories(epsilon INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Create version compatibility
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${PROJECT_BINARY_DIR}/epsilonConfigVersion.cmake"
    VERSION 0.1
    COMPATIBILITY AnyNewerVersion
)

# Install the library to the chosen CMAKE_INSTALL_PREFIX with normal folder structure
install(TARGETS epsilon
    EXPORT epsilonTargets
    LIBRARY DESTINATION lib COMPONENT Runtime
    ARCHIVE DESTINATION lib COMPONENT Development
    RUNTIME DESTINATION bin COMPONENT Runtime
    PUBLIC_HEADER DESTINATION include COMPONENT Development
    BUNDLE DESTINATION bin COMPONENT Runtime
)

# Create the configuration file to allow find_package(epsilon)
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/epsilonConfig.cmake.in"
    "${PROJECT_BINARY_DIR}/epsilonConfig.cmake"
    INSTALL_DESTINATION lib/cmake/mylib
)
install(EXPORT epsilonTargets DESTINATION lib/cmake/epsilon)
install(FILES "${PROJECT_BINARY_DIR}/epsilonConfigVersion.cmake"
              "${PROJECT_BINARY_DIR}/epsilonConfig.cmake"
        DESTINATION lib/cmake/epsilon)
# Copy over the header files
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION include)

if(EPSILON_ENABLE_UNIT_TESTS)
	message(STATUS "Building unit tests")
	# Unit test (does not get installed of course)
	add_executable(epsilon_unittest
		${PROJECT_SOURCE_DIR}/unittest/conversions.cpp
		${PROJECT_SOURCE_DIR}/unittest/eiconfig.hpp
		${PROJECT_SOURCE_DIR}/unittest/elementarytypes.cpp
		${PROJECT_SOURCE_DIR}/unittest/intersection2d.cpp
		${PROJECT_SOURCE_DIR}/unittest/intersection3d.cpp
		${PROJECT_SOURCE_DIR}/unittest/main.cpp
		${PROJECT_SOURCE_DIR}/unittest/matrix.cpp
		${PROJECT_SOURCE_DIR}/unittest/performance3d.hpp
		${PROJECT_SOURCE_DIR}/unittest/prime.cpp
		${PROJECT_SOURCE_DIR}/unittest/quaternion.cpp
		${PROJECT_SOURCE_DIR}/unittest/stdextensions.cpp
		${PROJECT_SOURCE_DIR}/unittest/types2d.cpp
		${PROJECT_SOURCE_DIR}/unittest/types3d.cpp
		${PROJECT_SOURCE_DIR}/unittest/unittest.cpp
		${PROJECT_SOURCE_DIR}/unittest/unittest.hpp
	)
	# The unit test makes some assumptions about include paths (not specifying 'ei/')
	target_include_directories(epsilon_unittest PRIVATE
		${PROJECT_SOURCE_DIR}/unittest
		${PROJECT_SOURCE_DIR}/include/ei
	)
	target_link_libraries(epsilon_unittest epsilon)
endif()