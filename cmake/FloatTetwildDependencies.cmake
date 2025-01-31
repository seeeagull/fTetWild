################################################################################
# Prepare dependencies
################################################################################

# Download external dependencies
include(FloatTetwildDownloadExternal)

################################################################################
# Required libraries
################################################################################

# fmt
if(NOT TARGET fmt::fmt)
	find_package(fmt REQUIRED)
endif()

# spdlog
if(NOT TARGET spdlog::spdlog)
	find_package(spdlog REQUIRED)

	# Create interface target
	add_library(spdlog INTERFACE)
	add_library(spdlog::spdlog ALIAS spdlog)
	target_link_libraries(spdlog INTERFACE fmt::fmt)
	target_compile_definitions(spdlog INTERFACE SPDLOG_FMT_EXTERNAL)
endif()

# Libigl
if(NOT TARGET igl::core)
	message(FATAL_ERROR "No libigl found. You probably need -DZENO_WITH_cgmesh:BOOL=ON")
endif()

# Geogram
if(NOT TARGET geogram::geogram)

	set(GEOGRAM_SEARCH_PATHS ${FLOAT_TETWILD_EXTERNAL}/geogram)
	find_path(GEOGRAM_SOURCE_INCLUDE_DIR
			geogram/basic/common.h
			PATHS ${GEOGRAM_SEARCH_PATHS}
			PATH_SUFFIXES src/lib
			NO_DEFAULT_PATH
	)
	set(GEOGRAM_ROOT ${GEOGRAM_SOURCE_INCLUDE_DIR}/../..)
	message(STATUS "Found Geogram here: ${GEOGRAM_ROOT}")

	if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		set(VORPALINE_ARCH_64 TRUE CACHE BOOL "" FORCE)
		set(VORPALINE_PLATFORM Win-vs-generic CACHE STRING "" FORCE)
		set(VORPALINE_BUILD_DYNAMIC false CACHE STRING "" FORCE)
	elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
		set(VORPALINE_PLATFORM Linux64-gcc CACHE STRING "" FORCE)
		set(VORPALINE_BUILD_DYNAMIC false CACHE STRING "" FORCE)
	elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		set(VORPALINE_PLATFORM Darwin-clang CACHE STRING "" FORCE)
		set(VORPALINE_BUILD_DYNAMIC false CACHE STRING "" FORCE)
	endif()

	option(GEOGRAM_WITH_GRAPHICS "Viewers and geogram_gfx library" OFF)
	option(GEOGRAM_WITH_LEGACY_NUMERICS "Legacy numerical libraries" OFF)
	option(GEOGRAM_WITH_HLBFGS "Non-linear solver (Yang Liu's HLBFGS)" OFF)
	option(GEOGRAM_WITH_TETGEN "Tetrahedral mesher (Hang Si's TetGen)" OFF)
	option(GEOGRAM_WITH_TRIANGLE "Triangle mesher (Jonathan Shewchuk's triangle)" OFF)
	option(GEOGRAM_WITH_EXPLORAGRAM "Experimental code (hexahedral meshing vpipeline and optimal transport)" OFF)
	option(GEOGRAM_WITH_LUA "Built-in LUA interpreter" OFF)
	option(GEOGRAM_LIB_ONLY "Libraries only (no example programs/no viewer)" ON)
	option(GEOGRAM_WITH_FPG "Predicate generator (Sylvain Pion's FPG)" OFF)
	option(GEOGRAM_USE_SYSTEM_GLFW3 "Use the version of GLFW3 installed in the system if found" OFF)

	add_subdirectory(${GEOGRAM_ROOT} geogram)
	target_include_directories(geogram SYSTEM PUBLIC ${GEOGRAM_SOURCE_INCLUDE_DIR})

	if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
		set_target_properties(geogram PROPERTIES COMPILE_FLAGS -fopenmp LINK_FLAGS -fopenmp)
		target_compile_options(geogram PUBLIC -fopenmp)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
		find_package(OpenMP REQUIRED)
		if(NOT TARGET OpenMP::OpenMP_CXX)
			add_library(OpenMP_TARGET INTERFACE)
				add_library(OpenMP::OpenMP_CXX ALIAS OpenMP_TARGET)
				target_compile_options(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
				find_package(Threads REQUIRED)
				target_link_libraries(OpenMP_TARGET INTERFACE Threads::Threads)
				target_link_libraries(OpenMP_TARGET INTERFACE ${OpenMP_CXX_FLAGS})
		endif()
		add_library(geogram_wrapper INTERFACE)
		add_library(geogram::geogram ALIAS geogram_wrapper)
		target_link_libraries(geogram_wrapper INTERFACE geogram OpenMP::OpenMP_CXX)
	else()
		add_library(geogram::geogram ALIAS geogram)
	endif()

	if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
		# remove warning for multiply defined symbols (caused by multiple
		# instantiations of STL templates)
		#target_compile_options(geogram INTERFACE /wd4251)

		# remove all unused stuff from windows.h
		# target_compile_definitions(geogram INTERFACE -DWIN32_LEAN_AND_MEAN)
		# target_compile_definitions(geogram INTERFACE -DVC_EXTRALEAN)

		# do not define a min() and a max() macro, breaks
		# std::min() and std::max() !!
		target_compile_definitions(geogram INTERFACE -DNOMINMAX)

		# we want M_PI etc...
		target_compile_definitions(geogram INTERFACE -D_USE_MATH_DEFINES)
	endif()

	add_library(geogram::geogram ALIAS geogram)
endif()

# C++11 threads
find_package(Threads REQUIRED)

# Json
if(NOT TARGET json)
	add_library(json INTERFACE)
	target_include_directories(json SYSTEM INTERFACE ${FLOAT_TETWILD_EXTERNAL}/json/include)
endif()

if(FLOAT_TETWILD_WITH_EXACT_ENVELOPE)
	float_tetwild_download_exact_envelope()
	add_subdirectory(${FLOAT_TETWILD_EXTERNAL}/exact_envelope)
endif()