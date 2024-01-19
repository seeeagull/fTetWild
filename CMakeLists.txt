################################################################################
cmake_minimum_required(VERSION 3.8)
project(FloatTetwild)
################################################################################

set(FLOAT_TETWILD_EXTERNAL ${CMAKE_CURRENT_SOURCE_DIR}/ext)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Color output
include(UseColors)

# Prepend function
include(PrependCurrentPath)

# FloatTetwild utils
include(FloatTetwildUtils)

# Extra warnings
include(Warnings)

# Sort projects inside the solution
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# Generate position independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################

# FloatTetwild options
option(FLOAT_TETWILD_ENABLE_TBB         "Enable TBB"                    ON)
option(FLOAT_TETWILD_USE_FLOAT          "Use floats instead of double"  ON)
option(FLOAT_TETWILD_WITH_EXACT_ENVELOPE    "Use exact envelope"        OFF)

# Sanitizer options
option(SANITIZE_ADDRESS   "Sanitize Address"       OFF)
option(SANITIZE_MEMORY    "Sanitize Memory"        OFF)
option(SANITIZE_THREAD    "Sanitize Thread"        OFF)
option(SANITIZE_UNDEFINED "Sanitize Undefined"     OFF)

# Options for libigl modules
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    OFF)
option(LIBIGL_WITH_PNG               "Use PNG"                      OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         ON)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)

# Setup dependencies
include(FloatTetwildDependencies)

################################################################################
# FloatTetwild library
################################################################################

# FloatTetWild API
add_library(FTetwildApi SHARED
	src/AABBWrapper.h src/AABBWrapper.cpp
	src/auto_table.hpp src/auto_table.cpp
	src/CSGTreeParser.hpp src/CSGTreeParser.cpp
	src/LocalOperations.h src/LocalOperations.cpp
	src/EdgeSplitting.h src/EdgeSplitting.cpp
	src/EdgeCollapsing.h src/EdgeCollapsing.cpp
	src/EdgeSwapping.h src/EdgeSwapping.cpp
	src/VertexSmoothing.h src/VertexSmoothing.cpp
	src/MeshImprovement.h src/MeshImprovement.cpp
	src/FloatTetwild.h src/FloatTetwild.cpp
	src/FloatTetDelaunay.h src/FloatTetDelaunay.cpp
	src/TriangleInsertion.cpp src/TriangleInsertion.h
	src/CutMesh.cpp src/CutMesh.h
	src/intersections.h src/intersections.cpp
	src/Simplification.h src/Simplification.cpp
	src/Mesh.cpp src/Mesh.hpp
	src/MeshIO.hpp src/MeshIO.cpp
	src/Logger.hpp src/Logger.cpp
	src/Parameters.h
	src/Types.hpp
	src/Statistics.h

	src/external/Predicates.cpp src/external/Predicates.hpp
	src/external/mesh_AABB.h src/external/mesh_AABB.cpp
	src/external/MshSaver.h src/external/MshSaver.cpp
	src/external/MshLoader.h src/external/MshLoader.cpp
	src/external/Exception.h
	src/external/triangle_triangle_intersection.cpp
	src/external/getRSS.c
	src/external/get_mem.h src/external/get_mem.cpp
	src/external/bfs_orient.h src/external/bfs_orient.cpp
	src/external/Rational.h
	src/external/WindingNumber.h
	
	api.h api.cpp)
target_compile_definitions(FTetwildApi PRIVATE -DFTetWild_EXPORT)
target_link_libraries(FTetwildApi PRIVATE tbb_from_instant_meshes)

find_package(GMPfTetWild)
IF(NOT ${GMP_FOUND})
	MESSAGE(FATAL_ERROR "Cannot find GMP")
ENDIF()

# find_package(MPFR)
# IF(NOT ${MPFR_FOUND})
# 	MESSAGE(FATAL_ERROR "Cannot find MPFR")
# ENDIF()

# add_library() can only be called without any source since CMake 3.11 ...
add_library(${PROJECT_NAME} src/Logger.cpp src/external/WindingNumber.h)

# Public include directory for FloatTetwild
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_BINARY_DIR}/include)
target_include_directories(FTetwildApi PUBLIC ${PROJECT_BINARY_DIR}/include)

# Extra warnings
target_link_libraries(${PROJECT_NAME} PRIVATE warnings::all)
target_link_libraries(FTetwildApi PRIVATE warnings::all)

# Use C++11
target_compile_features(${PROJECT_NAME} PUBLIC cxx_std_11)
target_compile_features(FTetwildApi PUBLIC cxx_std_11)
# if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
# 	target_compile_options(${PROJECT_NAME} PUBLIC "/Zc:__cplusplus")
# endif()
# target_compile_definitions(${PROJECT_NAME} PUBLIC -DEIGEN_STACK_ALLOCATION_LIMIT=8388608)

################################################################################
# Required libraries
################################################################################

if(FLOAT_TETWILD_USE_FLOAT)
	target_compile_definitions(${PROJECT_NAME} PUBLIC -DFLOAT_TETWILD_USE_FLOAT)
	target_compile_definitions(FTetwildApi PUBLIC -DFLOAT_TETWILD_USE_FLOAT)
endif()

target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${GMP_INCLUDE_DIRS})
target_include_directories(FTetwildApi SYSTEM PUBLIC ${GMP_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME}
	PUBLIC
		igl::core
		igl::predicates
		geogram::geogram
		spdlog::spdlog
		Threads::Threads
		json
		${GMP_LIBRARIES}
		# ${MPFR_LIBRARIES}
)
target_link_libraries(FTetwildApi
	PUBLIC
		igl::core
		igl::predicates
		geogram::geogram
		spdlog::spdlog
		Threads::Threads
		json
		${GMP_LIBRARIES}
)
include_directories(
  # geogram
  ${FLOAT_TETWILD_EXTERNAL}/geogram/src/lib
  # json
  ${FLOAT_TETWILD_EXTERNAL}/json/include
)
if(FLOAT_TETWILD_ENABLE_TBB)
	# find_package(TBB REQUIRED)
	target_link_libraries(${PROJECT_NAME} PRIVATE tbb_from_instant_meshes)
	target_compile_definitions(${PROJECT_NAME} PUBLIC FLOAT_TETWILD_USE_TBB)
	target_link_libraries(FTetwildApi PRIVATE tbb_from_instant_meshes)
	target_compile_definitions(FTetwildApi PUBLIC FLOAT_TETWILD_USE_TBB)
endif()

if(FLOAT_TETWILD_WITH_EXACT_ENVELOPE)
	target_link_libraries(${PROJECT_NAME} PUBLIC FastEnvelope)
	target_compile_definitions(${PROJECT_NAME} PUBLIC NEW_ENVELOPE)
	target_link_libraries(FTetwildApi PUBLIC FastEnvelope)
	target_compile_definitions(FTetwildApi PUBLIC NEW_ENVELOPE)
endif()

################################################################################
# Subdirectories
################################################################################

# Other sources
add_subdirectory(src)

