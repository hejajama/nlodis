configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gitsha1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp" @ONLY)

add_library(
	nlodisframework
    STATIC	
	nlodis.cpp
    nlodishelper_longitudinal.cpp
    nlodishelper_transverse.cpp
	qcd.cpp
    integration.cpp
    datatypes.cpp
    lo.cpp
	${CMAKE_CURRENT_BINARY_DIR}/gitsha1.cpp
)

add_executable(
    nlodis
    main.cpp
)

add_library(
    bkdipole
    STATIC
    dipole/bkdipole/bkdipole.cpp
    dipole/bkdipole/datafile.cpp
    dipole/bkdipole/interpolation.cpp
    dipole/dipoleamplitude.cpp
    dipole/vector.cpp
)


target_link_libraries(nlodisframework PUBLIC GSL::gsl GSL::gslcblas)
target_link_libraries(bkdipole PUBLIC GSL::gsl GSL::gslcblas)




add_executable(
    test
    tests/tests.cpp
    tests/test_massivedishelpers.cpp
    tests/test_nlodis_helpers.cpp
    tests/vector_tests.cpp
    )
	
include_directories(.)


# Build bundled Cuba library from vendored sources when requested
set(CUBA_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/Cuba-4.2.2")
set(CUBA_SRC_DIR "${CUBA_ROOT}/src")

# Build each algorithm in its own object library with precise include order
add_library(cuba_vegas_objs OBJECT ${CUBA_SRC_DIR}/vegas/Vegas.c)
target_include_directories(cuba_vegas_objs PRIVATE
    "${CUBA_ROOT}"
    "${CUBA_SRC_DIR}/vegas"
    "${CUBA_SRC_DIR}/common"
)
target_compile_definitions(cuba_vegas_objs PRIVATE NOUNDERSCORE=1)
set_target_properties(cuba_vegas_objs PROPERTIES C_STANDARD 99)

add_library(cuba_suave_objs OBJECT ${CUBA_SRC_DIR}/suave/Suave.c)
target_include_directories(cuba_suave_objs PRIVATE
    "${CUBA_ROOT}"
    "${CUBA_SRC_DIR}/suave"
    "${CUBA_SRC_DIR}/common"
)
target_compile_definitions(cuba_suave_objs PRIVATE NOUNDERSCORE=1)
set_target_properties(cuba_suave_objs PROPERTIES C_STANDARD 99)

add_library(cuba_divonne_objs OBJECT ${CUBA_SRC_DIR}/divonne/Divonne.c)
target_include_directories(cuba_divonne_objs PRIVATE
    "${CUBA_ROOT}"
    "${CUBA_SRC_DIR}/divonne"
    "${CUBA_SRC_DIR}/common"
)
target_compile_definitions(cuba_divonne_objs PRIVATE NOUNDERSCORE=1)
set_target_properties(cuba_divonne_objs PROPERTIES C_STANDARD 99)

add_library(cuba_cuhre_objs OBJECT ${CUBA_SRC_DIR}/cuhre/Cuhre.c)
target_include_directories(cuba_cuhre_objs PRIVATE
    "${CUBA_ROOT}"
    "${CUBA_SRC_DIR}/cuhre"
    "${CUBA_SRC_DIR}/common"
)
target_compile_definitions(cuba_cuhre_objs PRIVATE NOUNDERSCORE=1)
set_target_properties(cuba_cuhre_objs PROPERTIES C_STANDARD 99)

# Bundle into a single static library
add_library(cuba STATIC
    $<TARGET_OBJECTS:cuba_vegas_objs>
    $<TARGET_OBJECTS:cuba_suave_objs>
    $<TARGET_OBJECTS:cuba_divonne_objs>
    $<TARGET_OBJECTS:cuba_cuhre_objs>
    ${CUBA_SRC_DIR}/common/Data.c
    ${CUBA_SRC_DIR}/common/Global.c
)
target_include_directories(cuba PUBLIC "${CUBA_ROOT}")
target_compile_definitions(cuba PRIVATE NOUNDERSCORE=1)


set_target_properties(cuba PROPERTIES
    LINKER_LANGUAGE C
)


target_link_libraries(
	nlodis
	PRIVATE
    nlodisframework
    bkdipole
	cuba
	GSL::gsl
	GSL::gslcblas
)

target_link_libraries(
    test
    PRIVATE
    nlodisframework
    bkdipole
	cuba
	GSL::gsl
	GSL::gslcblas
)
