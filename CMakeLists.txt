# ###############################################################################
cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
set(CMAKE_CUDA_ARCHITECTURES 61 75)
project(MATTOPO LANGUAGES CXX CUDA)

# ###############################################################################
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/extern")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_DEBUG_POSTFIX "d")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")

# ###############################################################################
find_package(CUDA REQUIRED)
link_directories(${CUDA_TOOLKIT_ROOT_DIR}/lib/x64)

# find_package(OpenMP)
# if(OPENMP_FOUND)
#     set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

find_package(CUDA REQUIRED)
include_directories("${CUDA_INCLUDE_DIRS}")

find_package(ZLIB REQUIRED)

find_package(CGAL)

# Get Eigen3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
if (Eigen3_FOUND)
    # message(STATUS "${EIGEN3_VERSION_STRING}")
    include_directories(${EIGEN3_INCLUDE_DIR})
endif ()

# add dependencies
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
option(AUTO_DOWNLOAD "Auto download dependencies" ON)
include(mattopo_dependencies)

# ###############################################################################
# For MATTOPO binary
# ###############################################################################
set(HEADER_LIST
    "src/io_api.h"
)

set(SOURCE_LIST
    "src/main_voronoi.cpp" # yes, not .cu (for compiling CGAL)
    "src/io_api.cxx"
    "src/main_gui.cxx"
)

add_executable(${PROJECT_NAME} ${SOURCE_LIST} ${HEADER_LIST})
target_compile_options(${PROJECT_NAME} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    --generate-line-info
    --use_fast_math
    >)
# target_compile_features(${PROJECT_NAME} PUBLIC cuda_std_17)
target_compile_definitions(${PROJECT_NAME} PUBLIC -Dgeogram_EXPORTS)

# add common_headers include
target_include_directories(${PROJECT_NAME} PUBLIC
	$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>)

if(${MSVC})
    if(${CMAKE_VERSION} VERSION_GREATER 3.6.0)
        set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT "${PROJECT_NAME}")
    endif()
endif()

target_link_directories(${PROJECT_NAME} PUBLIC "${CMAKE_BINARY_DIR}/lib")
target_link_libraries(${PROJECT_NAME}
    general cuda
    general cublas
    optimized inputs
    optimized rpd3d
    optimized rpd3d_api
    optimized dist2mat
    optimized matbase
    optimized matfun
    optimized matfun_fix
    optimized IO
    optimized cpu_fallback # keep it
    optimized mshloader
    optimized polyscope
    optimized geogram
    optimized nlohmann_json::nlohmann_json
    optimized ZLIB::ZLIB
    optimized CGAL::CGAL
    optimized Eigen3::Eigen
    optimized igl::core
)
target_link_libraries(${PROJECT_NAME} general ${CUDA_LIBRARIES})

# add common headers
target_include_directories(${PROJECT_NAME} PUBLIC 
    ${EXTERNAL_DIR}/mat_modules/include
    ${EXTERNAL_DIR}/mat_modules/src/IO
    ${EXTERNAL_DIR}/mat_modules/src/matfun
    ${EXTERNAL_DIR}/mat_modules/src/matfun_fix
    ${EXTERNAL_DIR}/mat_modules/src/rpd3d_api
)