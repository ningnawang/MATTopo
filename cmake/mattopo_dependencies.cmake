# ###############################################################################
# CMake download helpers
# ###############################################################################

# download external dependencies
include(mattopo_downloadExternal)

# ###############################################################################
# Required dependencies
# ###############################################################################

if(NOT TARGET MAT_MODULES)
    rpd_download_mat_modules()
    add_subdirectory(${EXTERNAL_DIR}/mat_modules)
endif()