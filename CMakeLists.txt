cmake_minimum_required(VERSION 3.17)
project(amal
    VERSION 0.0.1
    LANGUAGES CXX C
    DESCRIPTION "App3d Math Library")

add_subdirectory(scripts)

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(${PROJECT_NAME}
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/public>
    $<INSTALL_INTERFACE:include>
)

# Internal
add_library(amal_internal INTERFACE)
add_dependencies(amal_internal amal_generated_headers)
target_include_directories(amal_internal INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/scripts
    )
target_link_libraries(${PROJECT_NAME} INTERFACE amal_internal)

set_target_properties(${PROJECT_NAME}
    PROPERTIES
    CXX_STANDARD 23
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS YES
)