#pragma once

#include <assert.h>
#include <geogram/mesh/mesh.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "common_geogram.h"
#include "io.h"
#include "io_utils.hpp"
#include "medial_mesh.h"
#include "medial_sphere.h"
#include "params.h"
#include "voronoi_defs.h"

void load_spheres_from_file(const char* filename,
                            std::vector<MedialSphere>& all_medial_spheres,
                            bool is_load_type, bool is_load_deleted,
                            bool is_save_dup_cnt);

void save_spheres_file(const std::vector<MedialSphere>& all_medial_spheres,
                       const std::string filename, bool is_save_type,
                       bool is_load_deleted, bool is_save_dup_cnt);

// // save ma
void export_ma_clean(const std::string& maname, const MedialMesh& mat);
void export_ma_ply(const std::string& maname, const MedialMesh& mat);

// save ma with dup_cnt and face importance
void export_ma_dup_cnt(const std::string& maname, MedialMesh& mat);
void load_ma_dup_cnt(const std::string& ma_path,
                     std::vector<MedialSphere>& all_medial_spheres,
                     MedialMesh& mat);

bool save_input_tet_houdini(const TetMesh& tet_mesh, std::string tet_name);