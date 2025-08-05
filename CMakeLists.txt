cmake_minimum_required(VERSION 3.17)
project(amal
    VERSION 0.0.1
    LANGUAGES CXX C
    DESCRIPTION "App3d Math Library")

add_subdirectory(scripts)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME}
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/scripts>
    $<INSTALL_INTERFACE:include>
)

add_dependencies(${PROJECT_NAME} amal_generated_headers)

set_target_properties(${PROJECT_NAME}
    PROPERTIES
    CXX_STANDARD 23
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS YES
)