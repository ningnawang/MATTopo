#ifndef H_MAIN_GUI_H
#define H_MAIN_GUI_H

#include <geogram/mesh/mesh.h>

#include <chrono>
#include <cmath>

#include "input_types.h"
#include "medial_mesh.h"
#include "medial_sphere.h"
#include "rpd_api.h"
#include "shrinking.h"
#include "triangulation.h"
#include "voronoi_defs.h"

class MainGuiWindow {
 private:
  static MainGuiWindow* instance_;
  static constexpr double point_radius_rel = 0.002;
  static constexpr double edge_radius_rel = 0.0008;
  using timer = std::chrono::system_clock;

 private:
  // non-pointers
  bool site_is_transposed;
  std::vector<float> site;
  std::vector<float> site_weights;
  std::vector<uint> site_flags;
  std::vector<int> site_knn;
  std::vector<float> site_cell_vol;
  std::vector<Vector3> pc_vertices;  // (non-restricted) powercell vertices
  std::set<aint2> checked_mmesh_edges;
  bool is_non_cad = false;
  // smaller the value, less aggresive the thinning
  double given_thinning_thres = 0.1;

  // pointers
  Parameter* params = nullptr;  // not const, update site_k
  // eg. const int * == int const * (pointer to const int)
  TetMesh* tet_mesh = nullptr;  // not const, init_corner_spheres() will update
  const SurfaceMesh* sf_mesh = nullptr;
  const std::set<int>* spheres_to_fix = nullptr;
  std::vector<MedialSphere>* all_medial_spheres = nullptr;  // not const
  //   std::vector<ConvexCellHost>* powercells = nullptr;     // not const
  MedialMesh* mmesh = nullptr;           // not const
  RegularTriangulationNN* rt = nullptr;  // not const
  RPD3D_GPU* rpd3d = nullptr;            // not const

  // For GEO FIX
  std::vector<double> fix_geo_samples;
  std::vector<float> fix_geo_samples_dist2mat;
  std::vector<aint3> fix_geo_samples_clostprim;

 public:
  void set_params(Parameter& _params);
  void set_tet_mesh(TetMesh& _tet_mesh);
  void set_tet_adj_info(const std::map<int, std::set<int>>& _v2tets,
                        const std::map<int, std::set<int>>& _tet_vs2sf_fids,
                        const std::map<aint4, int>& _tet_vs_lfs2tvs_map);
  void set_sf_mesh(const SurfaceMesh& _sf_mesh);
  void set_all_medial_spheres(std::vector<MedialSphere>& _all_medial_spheres);
  void set_rt(RegularTriangulationNN& _rt);
  void set_medial_mesh(MedialMesh& _mmesh);
  // call after set_tet_mesh(), set_sf_mesh(), set_params()
  void set_rpd3d(RPD3D_GPU& _rpd3d);
  void set_thinning_thres(const double& _thinning_thres);

 public:
  MainGuiWindow();
  ~MainGuiWindow();  // implemented in main_gui.cpp

  // show polyscope window
  void compute_partial_rpd_w_degeneration(
      std::vector<MedialSphere>& all_medial_spheres,
      const bool is_given_all_tets, bool is_debug);
  void compute_rpd_init();
  void show(bool is_compute_rpd = true);
  static void callbacks();

  // main iteration functions
  int num_itr_global = 0;
  int num_sphere_added = 0;

  // fix
  int num_topo_itr = 0;
  int run_topo_fix_itr(const SurfaceMesh& sf_mesh, const TetMesh& tet_mesh,
                       std::vector<MedialSphere>& all_medial_spheres,
                       std::vector<ConvexCellHost>& powercells,
                       bool is_debug = false);
  int num_geo_itr = 0;
  int run_fix_geo_itr(const SurfaceMesh& sf_mesh, const TetMesh& tet_mesh,
                      const Parameter& params, MedialMesh& mmesh,
                      std::vector<MedialSphere>& all_medial_spheres,
                      bool is_debug = false);
  int num_extf_itr = 0;
  int run_extf_fix_itr(const Parameter& param, TetMesh& tet_mesh,
                       const SurfaceMesh& sf_mesh,
                       const std::vector<ConvexCellHost>& powercells,
                       std::vector<MedialSphere>& all_medial_spheres,
                       MedialMesh& mmesh, bool is_debug = false);
  int num_intf_itr = 0;
  int run_intf_fix_itr(const SurfaceMesh& sf_mesh, const TetMesh& tet_mesh,
                       const MedialMesh& mmesh,
                       std::vector<MedialSphere>& all_medial_spheres,
                       bool is_debug = false);

  // generate medial mesh
  int num_mmesh_itr = 0;
  void run_mmesh_generator(const SurfaceMesh& sf_mesh, MedialMesh& mmesh,
                           std::vector<MedialSphere>& all_medial_spheres,
                           bool is_trace_mstruc = false,
                           bool is_update_tan_pls = true,
                           bool is_debug = false);

  // Convex cells
  int given_cell_id = -1;
  void show_one_cell(const int given_cell_id,
                     adouble3 ran_color = {{-1., -1., -1.}}, int suffix = -1);
  int given_site_id = -1;
  int site_cnt = 0;
  void show_one_site(const int given_site_id);
  void clear_cells();
  void show_one_site_boundary(const int given_site_id);
  void clear_site_boundary();
  void show_result_convex_cells(std::vector<ConvexCellHost>& powercells,
                                bool is_slice_plane = false);

  // Tetmesh
  void show_tet_mesh(const TetMesh& tet_mesh);
  int given_tet_id = -1;
  void show_one_tet(const TetMesh& tet_mesh, const int& given_tet_id);

  // Surface
  int given_sf_face_id = -1;
  void show_surface_mesh(const GEO::Mesh& sf_mesh, const int& given_sf_face_id);
  void convert_to_show_edges(const std::vector<adouble3>& old_pos,
                             const std::set<aint2>& old_edges,
                             std::string curve_name);

  int given_sphere_id = -1;
  void show_pin_points(const std::vector<MedialSphere>& all_medial_spheres,
                       const int given_sphere_id);
  void show_all_medial_spheres(
      const std::vector<MedialSphere>& all_medial_spheres);
  void show_one_sphere(const std::vector<MedialSphere>& all_medial_spheres,
                       const int given_sphere_id,
                       const bool is_show_radius = false,
                       const bool is_clear_all = false);

  // Medial Mesh
  int given_mat_face_id = -1;
  void show_medial_mesh(const MedialMesh& mmesh, int given_mat_face_id = -1);
  void show_medial_edges(const MedialMesh& mmesh);
  // debug
  int given_mat_edge_v1 = -1, given_mat_edge_v2 = -1;
  void show_medial_edge_info(MedialMesh& mmesh, int given_mat_edge_v1,
                             int given_mat_edge_v2);

  // Fix GEO
  int given_sample_id = -1;
  void show_fix_geo_samples(const std::vector<double>& geo_samples,
                            const std::vector<float>& samples_dist2mat,
                            const std::vector<aint3>& samples_clostprim,
                            const int given_sample_id);

  // Concave line pins
  void show_fe_samples(const std::vector<FeatureLine>& feature_lines);

  // automation
  void auto_all();
  void auto_all_non_cad();
  void auto_fix_topo();
  void auto_fix_extf();
  void auto_fix_intf();
  void auto_fix_geo();
};

#endif  // __H_MAIN_GUI_H__
