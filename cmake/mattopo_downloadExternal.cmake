# ###############################################################################
include(DownloadProject)
include(FetchContent)

# Shortcut function
function(rpd_download_project name)
	if(AUTO_DOWNLOAD)
		download_project(
			PROJ ${name}
			SOURCE_DIR ${EXTERNAL_DIR}/${name}
			DOWNLOAD_DIR ${EXTERNAL_DIR}/.cache/${name}
			${ARGN}
		)
	endif()
endfunction()

# ###############################################################################

# libigl
function(rpd_download_libigl)
	FetchContent_Declare(libigl
		GIT_REPOSITORY https://github.com/libigl/libigl.git
		GIT_TAG v2.5.0
	)
	FetchContent_MakeAvailable(libigl)
endfunction()

# geogram
function(tetwild_download_geogram)
	rpd_download_project(geogram
		GIT_REPOSITORY https://github.com/alicevision/geogram.git
		GIT_TAG v1.7.5
	)
endfunction()

## fmt
function(tetwild_download_fmt)
rpd_download_project(fmt
GIT_REPOSITORY https://github.com/fmtlib/fmt.git
GIT_TAG        5.2.0
)
endfunction()

## spdlog
function(tetwild_download_spdlog)
rpd_download_project(spdlog
GIT_REPOSITORY https://github.com/gabime/spdlog.git
GIT_TAG        v1.1.0
)
endfunction()

## CLI11
function(tetwild_download_cli11)
rpd_download_project(cli11
GIT_REPOSITORY     https://github.com/CLIUtils/CLI11
GIT_TAG            v1.6.1
)
endfunction()

# polyscope
function(rpd_download_polyscope)
	rpd_download_project(polyscope
		GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
		GIT_TAG eb07f8acaf5c8dba30de9587f84b14ba0c411ca1
	)
endfunction()

# json
function(rpd_download_json)
FetchContent_Declare(
	json 
	URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz
)
FetchContent_MakeAvailable(json)
endfunction()

# libmat
function(rpd_download_libmat)
	rpd_download_project(libmat
		GIT_REPOSITORY https://github.com/ningnawang/libmat.git
		GIT_TAG v0.0.1
	)
endfunction()