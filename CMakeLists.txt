cmake_minimum_required(VERSION 3.17)

if(NOT DEFINED ACBT_LOADED)
    include(cmake/utils.cmake)
endif()

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

option(AMAL_INTEGRATE_ACUL "Integrate acul" ON)
if(AMAL_INTEGRATE_ACUL)
    target_compile_definitions(${PROJECT_NAME} INTERFACE -DAMAL_INTEGRATE_ACUL)
    target_link_libraries(${PROJECT_NAME} INTERFACE acul)
endif()

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

set_target_properties(${PROJECT_NAME}
    PROPERTIES
    CXX_STANDARD 23
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS YES
)