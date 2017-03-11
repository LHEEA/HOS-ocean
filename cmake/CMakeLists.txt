CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
PROJECT(HOS_ocean)

SET(CMAKE_VERBOSE_MAKEFILE OFF)
SET(BUILD_SHARED_LIBS ON)

# Define manually these variables to access fftw and lapack libraries
# FFTW_ROOT
# LAPACK_ROOT

################################################################################
SET(SRC_DIRECTORY     ${CMAKE_CURRENT_SOURCE_DIR}/../sources)
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

IF(NOT DEFINED CMAKE_COLOR_MAKEFILE)
    SET(CMAKE_COLOR_MAKEFILE OFF)
ENDIF()

IF(CMAKE_SIZEOF_VOID_P EQUAL 8)
    SET(ARCH64_BITS TRUE)
ELSE()
    SET(ARCH64_BITS FALSE)
ENDIF()

IF(NOT DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Debug)
ENDIF()

MESSAGE(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
STRING(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
IF(CMAKE_BUILD_TYPE_UPPER MATCHES DEBUG)
    SET(BUILD_SHARED_LIBS OFF)
ENDIF()

IF(DEFINED INSTALL_PREFIX)
    SET(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX})
ELSE()
    SET(CMAKE_INSTALL_PREFIX ../install${CMAKE_BUILD_TYPE})
ENDIF()

################################################################################
IF(NOT DEFINED FFTW_USE_STATIC_LIBS)
    SET(FFTW_USE_STATIC_LIBS TRUE)
ENDIF()
FIND_PACKAGE(FFTW REQUIRED)
MESSAGE(STATUS "FFTW_INCLUDES  : ${FFTW_INCLUDES}")
MESSAGE(STATUS "FFTW_LIBRARIES : ${FFTW_LIBRARIES}")
MESSAGE(STATUS "FFTW_LIB       : ${FFTW_LIB}")
MESSAGE(STATUS "FFTWF_LIB      : ${FFTWF_LIB}")
MESSAGE(STATUS "FFTWL_LIB      : ${FFTWL_LIB}")

IF(DEFINED LAPACK_ROOT)
    FIND_PACKAGE(LAPACK REQUIRED
                 HINTS ${LAPACK_ROOT})
ELSE()
    FIND_PACKAGE(LAPACK REQUIRED)
ENDIF()
################################################################################

MESSAGE(STATUS "CMAKE_GENERATOR = ${CMAKE_GENERATOR}")
MESSAGE(STATUS "CMAKE_INSTALL_PREFIX = ${CMAKE_INSTALL_PREFIX}")

INCLUDE(CMakeOutputDirectories)

ENABLE_LANGUAGE(Fortran)
INCLUDE(CMakeCompilationOptionsFortran)

################################################################################
INCLUDE(CMakeSources.cmake)
################################################################################

MESSAGE(STATUS "CMAKE_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

ADD_EXECUTABLE(HOS_ocean ${srcHOS_ocean})
TARGET_LINK_LIBRARIES(HOS_ocean lapack blas ${FFTW_LIBRARIES})

ADD_EXECUTABLE(Post_processing ${srcPost_processing})
TARGET_LINK_LIBRARIES(Post_processing ${FFTW_LIBRARIES})

ADD_EXECUTABLE(check_W_DYue ${srcCheck_W_DYue})
TARGET_LINK_LIBRARIES(check_W_DYue lapack blas ${FFTW_LIBRARIES})

INSTALL(TARGETS HOS_ocean Post_processing check_W_DYue
        RUNTIME DESTINATION ${RUNTIME_OUTPUT_DIRECTORY}
        LIBRARY DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
SET(ADDITIONAL_LIB)
LIST(APPEND ADDITIONAL_LIB blas lapack)
FOREACH(f ${ADDITIONAL_LIB})
    IF(EXISTS ${f})
        INSTALL(FILES ${f} DESTINATION ${LIBRARY_OUTPUT_DIRECTORY})
    ENDIF()
ENDFOREACH()

ENABLE_TESTING()
FILE(WRITE "${CMAKE_CURRENT_BINARY_DIR}/input_HOS.dat"
"Choice of computed case      :: i_case           :: 3
--- Geometry of the horizontal domain
Length in x-direction        :: xlen             :: 20.0
Length in y-direction        :: ylen             :: 20.0
--- Time stuff
Duration of the simulation   :: T_stop           :: 1.0
Sampling frequency (output)  :: f_out            :: 1.0
Tolerance of RK scheme       :: toler            :: 1.0e-7
Dommermuth initialisation    :: n                :: 4
Dommermuth initialisation    :: Ta               :: 0.0
--- Physical dimensional parameters
Gravity                      :: grav             :: 9.81
Water depth                  :: depth            :: 35.
--- Irregular waves (i_case=3)
Peak period in s             :: Tp_real          :: 10.0
Significant wave height in m :: Hs_real          :: 4.5
JONSWAP Spectrum             :: gamma            :: 3.3
Directionality (Dysthe)      :: beta             :: 0.785398
--- Output files
Tecplot version              :: tecplot          :: 11
Output: 1-dim. ; 0-nondim.   :: i_out_dim        :: 1
3d free surface quantities   :: i_3d             :: 1
3d modes                     :: i_a_3d           :: 1
2d free surface, center line :: i_2d             :: 0
Wave probes in domain        :: i_prob           :: 0
Swense output 1='yes',0='no' :: i_sw             :: 0
")

FILE(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Results")
ADD_TEST(NAME first_test
         COMMAND HOS_ocean
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})


INCLUDE(CMakePack)
