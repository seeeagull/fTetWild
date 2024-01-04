################################################################################

include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(FLOAT_TETWILD_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(FLOAT_TETWILD_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(float_tetwild_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${FLOAT_TETWILD_EXTERNAL}/${name}
        DOWNLOAD_DIR ${FLOAT_TETWILD_EXTERNAL}/.cache/${name}
        QUIET
        ${FLOAT_TETWILD_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

## libigl
function(float_tetwild_download_libigl)
    float_tetwild_download_project(libigl
        GIT_REPOSITORY https://github.com/libigl/libigl.git
        GIT_TAG        abd9f0f365545652dbc39cea94f045d9dd99a10d
    )
endfunction()

## Json
function(float_tetwild_download_json)
    float_tetwild_download_project(json
        GIT_REPOSITORY https://github.com/jdumas/json
        GIT_TAG        0901d33bf6e7dfe6f70fd9d142c8f5c6695c6c5b
    )
endfunction()

## Geogram LGPL
function(float_tetwild_download_geogram)
    float_tetwild_download_project(geogram
        GIT_REPOSITORY https://github.com/seeeagull/geogram.git
        GIT_TAG        7705a929db874a2fdbe08d3a96e476c626b4ff53
    )
endfunction()

## exact envelope
function(float_tetwild_download_exact_envelope)
    float_tetwild_download_project(exact_envelope
            GIT_REPOSITORY https://github.com/wangbolun300/fast-envelope
            GIT_TAG        520ee04b6c69a802db31d1fd3a3e6e382d10ef98
            )
endfunction()