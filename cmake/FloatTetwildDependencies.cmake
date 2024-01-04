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
	# target_include_directories(spdlog INTERFACE ${FLOAT_TETWILD_EXTERNAL}/spdlog/include)
	target_link_libraries(spdlog INTERFACE fmt::fmt)
	target_compile_definitions(spdlog INTERFACE SPDLOG_FMT_EXTERNAL)
endif()

# Libigl
if(NOT TARGET igl::core)
	float_tetwild_download_libigl()

	# Import libigl targets
	list(APPEND CMAKE_MODULE_PATH "${FLOAT_TETWILD_EXTERNAL}/libigl/cmake")
	include(libigl)
endif()

# Geogram
if(NOT TARGET geogram::geogram)
	float_tetwild_download_geogram()
	include(geogram)
endif()


# TBB

# C++11 threads
find_package(Threads REQUIRED)


# Json
if(NOT TARGET json)
	float_tetwild_download_json()
	add_library(json INTERFACE)
	target_include_directories(json SYSTEM INTERFACE ${FLOAT_TETWILD_EXTERNAL}/json/include)
endif()

if(FLOAT_TETWILD_WITH_EXACT_ENVELOPE)
	float_tetwild_download_exact_envelope()
	add_subdirectory(${FLOAT_TETWILD_EXTERNAL}/exact_envelope)
endif()