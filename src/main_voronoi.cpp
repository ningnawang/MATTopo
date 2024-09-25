#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include "fix_topo.h"
#include "input_types.h"
#include "io.h"
#include "io_api.h"
#include "main_gui.h"
#include "medial_mesh.h"
#include "params.h"
#include "rpd_api.h"
#include "sharp_feature_detection.h"
#include "shrinking.h"
#include "thinning.h"
#include "triangulation.h"
#include "voronoi.h"

void printDevProp() {
  int devCount;  // Number of CUDA devices
  cudaError_t err = cudaGetDeviceCount(&devCount);
  if (err != cudaSuccess) {
    std::cerr << "Failed to initialize CUDA / failed to count CUDA devices "
                 "(error code << "
              << cudaGetErrorString(err) << ")! [file: " << __FILE__
              << ", line: " << __LINE__ << "]" << std::endl;
    exit(1);
  }

  printf("CUDA Device Query...\n");
  printf("There are %d CUDA devices.\n", devCount);

  // Iterate through devices
  for (uint i = 0; i < devCount; ++i) {
    // Get device properties
    printf("\nCUDA Device #%d\n", i);
    cudaDeviceProp devProp;
    cudaGetDeviceProperties(&devProp, i);
    printf("Major revision number:         %d\n", devProp.major);
    printf("Minor revision number:         %d\n", devProp.minor);
    printf("Name:                          %s\n", devProp.name);
    printf("Total global memory:           %lu\n", devProp.totalGlobalMem);
    printf("Total shared memory per block: %lu\n", devProp.sharedMemPerBlock);
    printf("Total registers per block:     %d\n", devProp.regsPerBlock);
    printf("Warp size:                     %d\n", devProp.warpSize);
    printf("Maximum memory pitch:          %lu\n", devProp.memPitch);
    printf("Maximum threads per block:     %d\n", devProp.maxThreadsPerBlock);
    for (uint i = 0; i < 3; ++i)
      printf("Maximum dimension %d of block:  %d\n", i,
             devProp.maxThreadsDim[i]);
    for (uint i = 0; i < 3; ++i)
      printf("Maximum dimension %d of grid:   %d\n", i, devProp.maxGridSize[i]);
    printf("Clock rate:                    %d\n", devProp.clockRate);
    printf("Total constant memory:         %lu\n", devProp.totalConstMem);
    printf("Texture alignment:             %lu\n", devProp.textureAlignment);
    printf("Concurrent copy and execution: %s\n",
           (devProp.deviceOverlap ? "Yes" : "No"));
    printf("Number of multiprocessors:     %d\n", devProp.multiProcessorCount);
    printf("Kernel execution timeout:      %s\n",
           (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
  }
}

int main(int argc, char** argv) {
  // printDevProp();
  std::string instruction =
      "./bin/MATTOPO <tet_mesh.msh> <num_spheres> <cad/non_cad> <optional: "
      "thinning_threshold, default=0.3, smaller less aggressive thinning>";
  if (4 > argc) {
    std::cerr << "Usage: " << argv[0] << instruction << std::endl;
    return 1;
  }

  std::srand(RAN_SEED);  // set random seed
  std::string tet_mesh_file = argv[1];
  uint init_num_spheres = atoi(argv[2]);

  // loading extra file
  // 0: no file
  // 1: sph file for spheres
  // 2: ma file for medial mesh
  // 3: non_cad model
  int load_file_flag = 0;
  std::string file_path;
  if (argv[3]) {
    file_path = argv[3];
    if (file_path == "non_cad") {
      printf("loading non_cad model: %s\n", tet_mesh_file.c_str());
      load_file_flag = 3;
    } else if (file_path == "cad") {
      std::string file_ext = get_file_ext(file_path);
      if (file_ext == "sph") {
        printf("loading sph file: %s\n", file_path.c_str());
        load_file_flag = 1;
      } else if (file_ext == "ma_dup") {
        printf("loading ma file: %s\n", file_path.c_str());
        load_file_flag = 2;
      }
    } else {
      std::cerr << "Usage: " << argv[0] << instruction << std::endl;
      return 1;
    }
  }

  int* initptr = nullptr;
  cudaError_t err = cudaMalloc(
      &initptr, sizeof(int));  // unused memory, needed for initialize the
                               // GPU before time measurements
  if (err != cudaSuccess) {
    std::cerr << "Failed to allocate (error code << " << cudaGetErrorString(err)
              << ")! [file: " << __FILE__ << ", line: " << __LINE__ << "]"
              << std::endl;
    return 1;
  }

  // Load tet from file
  TetMesh tet_mesh(tet_mesh_file);
  Parameter params;
  // please set normalize to true for model joint.tet/.xyz
  if (!load_tet(tet_mesh_file, tet_mesh.tet_vertices, tet_mesh.tet_indices,
                true /*normalize*/, params)) {
    std::cerr << tet_mesh_file << ": could not load file" << std::endl;
    return 1;
  }
  load_v2tets(tet_mesh.tet_vertices, tet_mesh.tet_indices, tet_mesh.v2tets);
  // for (int i = 0; i < tet_mesh.tet_indices.size() / 4; i++) {
  //   if (tet_mesh.tet_indices[i * 4] == 768 ||
  //       tet_mesh.tet_indices[i * 4 + 1] == 768 ||
  //       tet_mesh.tet_indices[i * 4 + 2] == 768 ||
  //       tet_mesh.tet_indices[i * 4 + 3] == 768) {
  //     printf("tet %d has tvids (%d,%d,%d,%d) \n", i,
  //            tet_mesh.tet_indices[i * 4], tet_mesh.tet_indices[i * 4 + 1],
  //            tet_mesh.tet_indices[i * 4 + 2], tet_mesh.tet_indices[i * 4 +
  //            3]);
  //   }
  // }

  // Read surface file .geogram if exists
  // we store the mapping from tet_vertices to sf_mesh in vertex attributes
  //
  // this will make sure all geogram algoriths can be used
  GEO::initialize();
  GEO::CmdLine::import_arg_group("algo");
  std::string surface_path = get_other_file_path(tet_mesh_file, 0 /*.geogram*/);
  std::cout << "surface path: " << surface_path << std::endl;
  SurfaceMesh sf_mesh;
  if (!is_file_exist(surface_path)) {
    std::cout << "surface path: " << surface_path
              << " not exist, load from tet_mesh_file: " << tet_mesh_file
              << std::endl;
    std::vector<std::array<float, 3>> sf_vertices;
    std::vector<std::array<int, 3>> sf_faces;
    get_surface_from_tet(tet_mesh.tet_vertices, tet_mesh.tet_indices,
                         sf_vertices, sf_faces, sf_mesh.sf2tet_vs_mapping);
    load_sf_mesh_from_internal(sf_vertices, sf_faces, sf_mesh.sf2tet_vs_mapping,
                               sf_mesh);
    if (!save_sf_mesh_geogram(surface_path, sf_mesh)) {
      std::cerr << surface_path << ": could not save surface file" << std::endl;
      return 1;
    }
    // for debug
    // save surface mesh as obj file
    save_sf_mesh(get_other_file_path(tet_mesh_file, 1 /*.obj*/), sf_mesh);
  } else {
    if (!load_surface_mesh_geogram(surface_path, sf_mesh)) {
      std::cerr << surface_path << ": could not load surface file" << std::endl;
      return 1;
    }
  }
  // // save scaled surface mesh as obj file
  // save_sf_mesh_scaled(get_other_file_path(tet_mesh_file, 3 /*.obj*/),
  // sf_mesh, params);

  std::cout << "loaded surface_mesh #v: " << sf_mesh.vertices.nb()
            << ", #f: " << sf_mesh.facets.nb() << std::endl;
  // ATTENTION: this will REORDER GEO::Mesh!!!!
  pre_and_init_aabb(sf_mesh, sf_mesh.aabb_wrapper);
  // load the mapping from tet_vertices to sf_mesh
  // NOTE: must use after pre_and_init_aabb(), which will REORDER GEO::Mesh
  std::map<int, int> _;
  std::map<int, std::set<int>> __;
  load_sf_tet_mapping(sf_mesh, _, __, tet_mesh.tet_vs2sf_fids);
  load_tet_adj_info(tet_mesh.v2tets, tet_mesh.tet_indices,
                    tet_mesh.tet_vs2sf_fids, sf_mesh.facets.nb(),
                    tet_mesh.v_adjs, tet_mesh.e_adjs, tet_mesh.f_adjs,
                    tet_mesh.f_ids);
  sf_mesh.reload_sf2tet_vs_mapping();
  sf_mesh.reload_sf_fid_neighs();

  // Detect sharp features
  // (store in tet mesh indices, matching tet_vertices)
  detect_mark_sharp_features(params, sf_mesh, tet_mesh);
  pre_and_init_feature_aabb(tet_mesh, sf_mesh.aabb_wrapper);
  sf_mesh.reload_sf_fid_neighs_no_cross();

  // save sf_mesh with extf marked in .ma format
  save_sf_mesh_with_extf(get_other_file_path(tet_mesh_file, 4 /*.ma*/),
                         sf_mesh);

  // Init medial spheres
  std::vector<MedialSphere> all_medial_spheres;
  if (load_file_flag == 1) {
    load_spheres_from_file(file_path.data(), all_medial_spheres,
                           true /*is_load_type*/, true /*is_load_deleted*/,
                           true /*is_save_dup_cnt*/);
  } else if (load_file_flag == 3) {  // non_cad
    sf_mesh.is_non_cad = true;
    tet_mesh.is_non_cad = true;
    init_and_shrink(sf_mesh, tet_mesh, all_medial_spheres, init_num_spheres,
                    -1 /*itr_limit*/, false /*is_debug*/);
  } else {
    init_and_shrink(sf_mesh, tet_mesh, all_medial_spheres, init_num_spheres,
                    -1 /*itr_limit*/, false /*is_debug*/);
    insert_spheres_for_concave_lines_new(
        sf_mesh, tet_mesh.cc_corners, tet_mesh.feature_edges, tet_mesh.ce_lines,
        all_medial_spheres,
        params.cc_len_eps_rel * params.bbox_diag_l /*cc_len_eps*/,
        false /*is_debug*/);
  }

  RegularTriangulationNN rt;
  RPD3D_GPU rpd3d;
  MedialMesh mmesh;
  if (load_file_flag == 1) mmesh.is_corner_spheres_created = true;

  // show Gui
  MainGuiWindow main_gui;
  main_gui.set_params(params);
  main_gui.set_tet_mesh(tet_mesh);
  main_gui.set_sf_mesh(sf_mesh);
  main_gui.set_all_medial_spheres(all_medial_spheres);
  main_gui.set_rt(rt);
  main_gui.set_medial_mesh(mmesh);
  main_gui.set_rpd3d(rpd3d);

  // set threshold for thinning
  main_gui.set_thinning_thres(0.1);
  if (argv[4]) {
    float thres = atof(argv[4]);
    if (thres > 0) {
      printf("Setting thinning threshold to: %f\n", thres);
      main_gui.set_thinning_thres(thres);
    }
  }

  if (load_file_flag == 1) {
    // load spheres .sph
    main_gui.compute_rpd_init();
    mmesh.validate_mmesh_dup_cnt();
    // prune(all_medial_spheres, mmesh, main_gui.given_thinning_thres,
    //       true /*is_dup_cnt*/, false /*is_import_given*/,
    //       false /*is_sort_randomly*/, false /*is_debug*/);
    // compute_Euler(mmesh);
    // std::string name_no_ext =
    //     get_only_file_name(tet_mesh.tet_path_with_ext, false
    //                        /*withExtension*/);
    // export_ma_ply(name_no_ext, mmesh);
    main_gui.show(false /*is_compute_rpd*/);
  } else if (load_file_flag == 2) {
    // load medial mesh .ma_dup
    load_ma_dup_cnt(file_path, all_medial_spheres, mmesh);
    mmesh.validate_mmesh_dup_cnt();
    main_gui.show(false /*is_compute_rpd*/);
  } else if (load_file_flag == 3) {
    printf("processing NON_CAD model\n");
    // for non_cad models
    params.thres_concave = 0.8;
    params.thres_convex = 40;
    main_gui.auto_all_non_cad();
    main_gui.show(false /*is_compute_rpd*/);
  } else {
    // for cad models
    printf("processing CAD model\n");
    main_gui.auto_all();
    main_gui.show(false /*is_compute_rpd*/);
  }

  cudaFree(initptr);

  return 0;
}
