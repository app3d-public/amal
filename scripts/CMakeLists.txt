find_package(Python3 REQUIRED COMPONENTS Interpreter)

set(MATRIX_DEF_TYPE_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/matrix_type.py")
set(MATRIX_TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/templates/matrix.tpp.j2")
set(VECTOR_DEF_TYPE_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/vector_type.py")
set(VECTOR_TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/templates/vector.tpp.j2")

set(GENERATED_HEADERS "")

foreach(C RANGE 2 4)
    set(VECTOR_OUT_FILE "${CMAKE_CURRENT_BINARY_DIR}/amal/internal/vec${C}.hpp")

    add_custom_command(
        OUTPUT ${VECTOR_OUT_FILE}
        COMMAND ${CMAKE_COMMAND} -E echo "Generating ${VECTOR_OUT_FILE}"
        COMMAND ${Python3_EXECUTABLE} ${VECTOR_DEF_TYPE_SCRIPT} ${C} > ${VECTOR_OUT_FILE}
        DEPENDS ${VECTOR_DEF_TYPE_SCRIPT} ${VECTOR_TEMPLATE}
        VERBATIM
        )

    list(APPEND GENERATED_HEADERS ${VECTOR_OUT_FILE} ${MATRIX_OUT_FILE})

    foreach(R RANGE 2 4)
        set(MATRIX_OUT_FILE "${CMAKE_CURRENT_BINARY_DIR}/amal/internal/mat${C}x${R}.hpp")

        add_custom_command(
            OUTPUT ${MATRIX_OUT_FILE}
            COMMAND ${CMAKE_COMMAND} -E echo "Generating ${MATRIX_OUT_FILE}"
            COMMAND ${Python3_EXECUTABLE} ${MATRIX_DEF_TYPE_SCRIPT} ${C} ${R} > ${MATRIX_OUT_FILE}
            DEPENDS ${MATRIX_DEF_TYPE_SCRIPT} ${MATRIX_TEMPLATE}
            VERBATIM
        )

        list(APPEND GENERATED_HEADERS ${MATRIX_OUT_FILE})
    endforeach()
endforeach()

set(MATRIX_MULTIPLY_TEMPLATE   "${CMAKE_CURRENT_SOURCE_DIR}/templates/matrix_multiply.tpp.j2")
set(MATRIX_MULTIPLY_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/matrix_multiply.py")
set(MATRIX_MULTIPLY_OUT "${CMAKE_CURRENT_BINARY_DIR}/amal/internal/matrix_multiply.hpp")

add_custom_command(
    OUTPUT ${MATRIX_MULTIPLY_OUT}
    COMMAND ${CMAKE_COMMAND} -E echo "Generating ${MATRIX_MULTIPLY_OUT}"
    COMMAND ${Python3_EXECUTABLE} ${MATRIX_MULTIPLY_SCRIPT} > ${MATRIX_MULTIPLY_OUT}
    DEPENDS ${MATRIX_MULTIPLY_SCRIPT} ${MATRIX_MULTIPLY_TEMPLATE}
    VERBATIM
)
list(APPEND GENERATED_HEADERS ${MATRIX_MULTIPLY_OUT})

add_custom_target(amal_generated_headers DEPENDS ${GENERATED_HEADERS})