# ###############################################################################
# CMake download helpers
# ###############################################################################

# download external dependencies
include(mattopo_downloadExternal)

# ###############################################################################
# Required dependencies
# ###############################################################################

if(NOT TARGET LIBMAT)
    rpd_download_libmat()
    add_subdirectory(${EXTERNAL_DIR}/libmat)
endif()