
#include "main_gui.h"

#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/volume_mesh.h>

#include "fix_extf.h"
#include "fix_geo.h"
#include "fix_geo_error.h"
#include "fix_intf.h"
#include "fix_topo.h"
#include "io.h"
#include "io_api.h"
#include "rpd_api.h"
#include "thinning.h"
#include "voronoi.h"

MainGuiWindow* MainGuiWindow::instance_ = nullptr;
MainGuiWindow::MainGuiWindow() {
  if (instance_ != nullptr) {
    printf("ERROR: GuiWindow instance is not nullptr!!");
    exit(1);
  }
  instance_ = this;
}

MainGuiWindow::~MainGuiWindow() {
  // we delete here since we allocate memory
  // (called new) for each of these variables
  delete params;
  delete tet_mesh;
  delete sf_mesh;
  delete all_medial_spheres;
  delete rpd3d;

  this->instance_ = nullptr;
}

void MainGuiWindow::set_params(Parameter& _params) { params = &_params; }
void MainGuiWindow::set_tet_mesh(TetMesh& _tet_mesh) {
  tet_mesh = &_tet_mesh;
  is_non_cad = tet_mesh->is_non_cad;
}
void MainGuiWindow::set_sf_mesh(const SurfaceMesh& _sf_mesh) {
  sf_mesh = &_sf_mesh;
  is_non_cad = sf_mesh->is_non_cad;
}
void MainGuiWindow::set_all_medial_spheres(
    std::vector<MedialSphere>& _all_medial_spheres) {
  all_medial_spheres = &_all_medial_spheres;
}
void MainGuiWindow::set_rt(RegularTriangulationNN& _rt) { rt = &_rt; }
void MainGuiWindow::set_medial_mesh(MedialMesh& _mmesh) { mmesh = &_mmesh; }
void MainGuiWindow::set_rpd3d(RPD3D_GPU& _rpd3d) {
  rpd3d = &_rpd3d;
  rpd3d->init(this->tet_mesh, this->sf_mesh, this->params);
  rpd3d->set_spheres(this->all_medial_spheres);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main Iterations functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int MainGuiWindow::run_topo_fix_itr(
    const SurfaceMesh& sf_mesh, const TetMesh& tet_mesh,
    std::vector<MedialSphere>& all_medial_spheres,
    std::vector<ConvexCellHost>& powercells, bool is_debug) {
  printf("--------------------Topo Check %d--------------------\n",
         this->num_topo_itr);

  std::set<int> spheres_to_fix;
  check_cc_and_euler(powercells, all_medial_spheres, spheres_to_fix,
                     is_debug /*is_debug*/);

  if (spheres_to_fix.empty()) {
    polyscope::info("Topology check has passed!");
    return 0;
  }

  // int old_size = all_medial_spheres.size();
  int num_sphere_change = fix_topo_by_adding_new_sphere(
      instance_->num_itr_global, this->num_topo_itr, powercells, sf_mesh,
      tet_mesh, spheres_to_fix, all_medial_spheres, is_debug /*is_debug*/);
  // int num_added = all_medial_spheres.size() - old_size;

  polyscope::info("Topo check change " + std::to_string(num_sphere_change) +
                  " spheres");
  this->num_topo_itr++;
  return num_sphere_change;
}

int MainGuiWindow::run_fix_geo_itr(
    const SurfaceMesh& sf_mesh, const TetMesh& tet_mesh,
    const Parameter& params, MedialMesh& mmesh,
    std::vector<MedialSphere>& all_medial_spheres, bool is_debug) {
  printf("--------------------Geo Check %d--------------------\n", num_geo_itr);

  if (mmesh.vertices == nullptr) {
    // generate mmesh first
    run_mmesh_generator(sf_mesh, mmesh, all_medial_spheres);
  }
  int old_size = all_medial_spheres.size();
  check_and_fix_mm_geo(instance_->num_itr_global, sf_mesh, tet_mesh, params,
                       mmesh, all_medial_spheres, this->fix_geo_samples,
                       this->fix_geo_samples_dist2mat,
                       this->fix_geo_samples_clostprim, is_debug /*is_debug*/);
  int num_added = all_medial_spheres.size() - old_size;

  polyscope::info("Geo check added " + std::to_string(num_added) + " spheres");
  num_geo_itr++;
  return num_added;
}

int MainGuiWindow::run_extf_fix_itr(
    const Parameter& param, TetMesh& tet_mesh, const SurfaceMesh& sf_mesh,
    const std::vector<ConvexCellHost>& powercells,
    std::vector<MedialSphere>& all_medial_spheres, MedialMesh& mmesh,
    bool is_debug) {
  printf("--------------------EXTF Check %d--------------------\n",
         num_extf_itr);

  // add conrers first
  if (!mmesh.is_corner_spheres_created && !tet_mesh.corners_se_tet.empty()) {
    int num_corners = init_corner_spheres(instance_->num_itr_global, tet_mesh,
                                          all_medial_spheres);
    mmesh.is_corner_spheres_created = true;
    polyscope::info("EXTF check add " + std::to_string(num_corners) +
                    " corners");
    num_extf_itr++;
    return num_corners;
  }

  int old_size = all_medial_spheres.size();
  int num_changed = check_and_fix_external_feature(
      instance_->num_itr_global, param, tet_mesh.tet_vertices, powercells,
      sf_mesh, tet_mesh.tet_vs_lfs2tvs_map, tet_mesh.fl2corner_sphere,
      all_medial_spheres, is_debug /*is_debug*/);
  int num_added = all_medial_spheres.size() - old_size;
  polyscope::info("EXTF check changed " + std::to_string(num_changed) +
                  " spheres");
  num_extf_itr++;
  return num_changed;
}

int MainGuiWindow::run_intf_fix_itr(
    const SurfaceMesh& sf_mesh, const TetMesh& tet_mesh,
    const MedialMesh& mmesh, std::vector<MedialSphere>& all_medial_spheres,
    bool is_debug) {
  printf("--------------------INTF Check %d--------------------\n",
         num_intf_itr);

  int old_size = all_medial_spheres.size();
  random_check_edges_and_fix_internal_feature(
      instance_->num_itr_global, sf_mesh, tet_mesh, mmesh, all_medial_spheres,
      this->checked_mmesh_edges, is_debug);
  int num_added = all_medial_spheres.size() - old_size;

  // update T2 by TN
  int num_updated = 0;
  polyscope::info("INTF check added " + std::to_string(num_added) +
                  ", updated " + std::to_string(num_updated) + " spheres");
  num_intf_itr++;
  return num_added + num_updated;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Generate medial mesh
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void MainGuiWindow::run_mmesh_generator(
    const SurfaceMesh& sf_mesh, MedialMesh& mmesh,
    std::vector<MedialSphere>& all_medial_spheres, bool is_trace_mstruc,
    bool is_update_tan_pls, bool is_debug) {
  printf("--------------------Cal MedialMesh--------------------\n");
  mmesh.clear();
  mmesh.genearte_medial_spheres(all_medial_spheres);
  mmesh.generate_medial_edges(sf_mesh, is_update_tan_pls);
  mmesh.generate_medial_faces();
  mmesh.compute_faces_st_meta(sf_mesh.aabb_wrapper);  // for fix_geo
  mmesh.generate_medial_tets();
  // if (is_trace_mstruc) mmesh.trace_medial_structure(true /*is_debug*/);

  // mmesh.check_and_store_unthin_tets_in_mat();
  // compute euler
  // int euler = mmesh.vertices->size() - mmesh.edges.size() +
  // mmesh.faces.size();
  int euler = compute_Euler(mmesh);
  printf("[Euler] MedialMesh has Euler %d: v %ld, e %ld, f %ld, t %ld \n",
         euler, mmesh.vertices->size(), mmesh.edges.size(), mmesh.faces.size(),
         mmesh.tets.size());
  polyscope::info("Done calculating medial mesh");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_convex_cell_tet(ConvexCellHost& cc_trans,
                         std::vector<cfloat3>& voro_points,
                         std::vector<std::array<int, 4>>& voro_tets,
                         std::vector<int>& voro_tets_sites,
                         std::vector<int>& voro_tets_cell_ids,
                         std::vector<float>& voro_tets_euler) {
  if (!cc_trans.is_pc_explicit) cc_trans.reload_pc_explicit();
  // assert(cc_trans.is_pc_explicit);
  if (cc_trans.is_vertex_null) {
    printf("[ERROR] do not show cell %d since some vertices are null \n",
           cc_trans.id);
    return;
  }

  // printf("saving cell of site %d \n", cc_trans.voro_id);
  int row = voro_points.size();  // index from 0
  int bary_id = -1;
  cfloat3 bary = {0.0f, 0.0f, 0.0f};  // barycenter of cell
  for (const auto& voro_vertex : cc_trans.pc_points) {
    bary = cplus3(bary, voro_vertex);
    voro_points.push_back(
        cmake_float3(voro_vertex.x, voro_vertex.y, voro_vertex.z));
  }
  bary = cdivide3(bary, cc_trans.nb_v);
  // voro_points.push_back(make_float3(bary.x, bary.y, bary.z));
  voro_points.push_back(bary);
  bary_id = voro_points.size() - 1;

  // some clipping planes may not exist in tri but
  // we still sotre it, here is to filter those planes
  if (!cc_trans.is_active_updated) cc_trans.reload_active();
  assert(cc_trans.is_active_updated);
  const std::vector<int>& active_clipping_planes =
      cc_trans.active_clipping_planes;

  std::vector<std::vector<int>> voro_local_faces;
  int lf = 0;  // local fid

  FOR(plane, cc_trans.nb_p) {
    if (active_clipping_planes[plane] > 0) {
      std::vector<int> tab_v;   // index of dual vertex
      std::vector<int> tab_lp;  // local index of dual vertex in dual triangle
      // for each dual triangle
      FOR(t, cc_trans.nb_v) {
        // store info of dual vertex
        if ((int)cc_trans.ver_trans(t).x == plane) {
          tab_v.push_back(t);
          tab_lp.push_back(0);
        } else if ((int)cc_trans.ver_trans(t).y == plane) {
          tab_v.push_back(t);
          tab_lp.push_back(1);
        } else if ((int)cc_trans.ver_trans(t).z == plane) {
          tab_v.push_back(t);
          tab_lp.push_back(2);
        }
      }

      if (tab_lp.size() <= 2) {
        std::cout << (int)plane << std::endl;
      }

      int i = 0;
      int j = 0;
      voro_local_faces.push_back(std::vector<int>(0));

      while (voro_local_faces[lf].size() < tab_lp.size()) {
        int ind_i = (tab_lp[i] + 1) % 3;
        bool temp = false;
        j = 0;
        while (temp == false) {
          int ind_j = (tab_lp[j] + 2) % 3;
          if ((int)cc_trans.ith_plane(tab_v[i], ind_i) ==
              (int)cc_trans.ith_plane(tab_v[j], ind_j)) {
            voro_local_faces[lf].push_back(tab_v[i]);
            temp = true;
            i = j;
          }
          j++;
        }
      }

      // triangulate face, then make a tet
      // (bary, p1, p2, p3) makes a tet
      int nb_pts = voro_local_faces[lf].size();
      for (uint p = 1; p < voro_local_faces[lf].size() - 1; p++) {
        // (0, p, p+1) is a triangulation of face
        voro_tets.push_back({{bary_id, row + voro_local_faces[lf][0],
                              row + voro_local_faces[lf][p],
                              row + voro_local_faces[lf][(p + 1) % nb_pts]}});
        voro_tets_sites.push_back(cc_trans.voro_id);
        voro_tets_euler.push_back(cc_trans.euler);
        voro_tets_cell_ids.push_back(cc_trans.id);
      }

      // voro_faces += "f";
      // FOR(i, result[ind].size()) {
      //   voro_faces += " ";
      //   voro_faces += std::to_string(row + voro_faces[ind][i]);
      // }
      // voro_faces += "\n";
      lf++;
    }
  }
}

// 3d triangle mesh
void get_convex_cell_faces(
    ConvexCellHost& cc_trans, std::vector<cfloat3>& voro_points,
    std::vector<float>& voro_points_num_adjs, std::vector<aint3>& voro_faces,
    std::vector<int>& voro_faces_sites, std::vector<int>& voro_faces_cell_ids,
    std::vector<int>& voro_faces_tet_ids,
    std::vector<float>& voro_faces_cell_eulers,
    std::vector<int>& voro_faces_ids, std::vector<int>& voro_faces_sf_ids,
    std::vector<int>& voro_faces_neigh_id, std::vector<int>& voro_faces_num_adj,
    int max_sf_fid, bool is_boundary_only) {
  // assert(cc_trans.is_pc_explicit);
  if (!cc_trans.is_pc_explicit) cc_trans.reload_pc_explicit();
  if (cc_trans.is_vertex_null) {
    printf("[ERROR] do not show cell %d since some vertices are null \n",
           cc_trans.id);
    return;
  }

  int row = voro_points.size();  // index from 0
  // save all vertices
  FOR(i, cc_trans.pc_points.size()) {
    voro_points.push_back(cc_trans.pc_points[i]);
    voro_points_num_adjs.push_back((int)cc_trans.ver_trans(i).w);
  }

  // some clipping planes may not exist in tri but
  // we still sotre it, here is to filter those planes
  if (!cc_trans.is_active_updated) cc_trans.reload_active();
  assert(cc_trans.is_active_updated);
  const std::vector<int>& active_clipping_planes =
      cc_trans.active_clipping_planes;

  FOR(plane, cc_trans.nb_p) {
    if (active_clipping_planes[plane] <= 0) continue;
    cint2 hp = cc_trans.clip_id2_const(plane);

    // check if only shown boundary or not
    // 1. halfplane; 2. on sf_mesh fid
    bool is_shown = true;
    if (is_boundary_only) {
      is_shown = false;
      if ((hp.y != -1) || (hp.y == -1 && hp.x < max_sf_fid)) {
        is_shown = true;
      }
    }
    if (!is_shown) continue;

    std::vector<int> tab_v;   // index of dual vertex
    std::vector<int> tab_lp;  // local index of dual vertex in dual triangle
    // for each dual triangle
    FOR(t, cc_trans.nb_v) {
      // store info of dual vertex
      if ((int)cc_trans.ver_trans(t).x == plane) {
        tab_v.push_back(t);
        tab_lp.push_back(0);
      } else if ((int)cc_trans.ver_trans(t).y == plane) {
        tab_v.push_back(t);
        tab_lp.push_back(1);
      } else if ((int)cc_trans.ver_trans(t).z == plane) {
        tab_v.push_back(t);
        tab_lp.push_back(2);
      }
    }

    if (tab_lp.size() <= 2) {
      std::cout << (int)plane << std::endl;
    }

    int i = 0;
    int j = 0;
    int lf = 0;  // local fid
    std::vector<std::vector<int>> voro_local_faces;
    voro_local_faces.push_back(std::vector<int>(0));

    while (voro_local_faces[lf].size() < tab_lp.size()) {
      int ind_i = (tab_lp[i] + 1) % 3;
      bool temp = false;
      j = 0;
      while (temp == false) {
        int ind_j = (tab_lp[j] + 2) % 3;
        if ((int)cc_trans.ith_plane(tab_v[i], ind_i) ==
            (int)cc_trans.ith_plane(tab_v[j], ind_j)) {
          voro_local_faces[lf].push_back(tab_v[i]);
          temp = true;
          i = j;
        }
        j++;
      }
    }

    // triangulate face, then make a tet
    // (bary, p1, p2, p3) makes a tet
    int nb_pts = voro_local_faces[lf].size();
    for (uint p = 1; p < voro_local_faces[lf].size() - 1; p++) {
      // (0, p, p+1) is a triangulation of face
      voro_faces.push_back(
          {{row + voro_local_faces[lf][0], row + voro_local_faces[lf][p],
            row + voro_local_faces[lf][(p + 1) % nb_pts]}});
      voro_faces_sites.push_back(cc_trans.voro_id);
      voro_faces_cell_eulers.push_back(cc_trans.euler);
      voro_faces_cell_ids.push_back(cc_trans.id);
      voro_faces_tet_ids.push_back(cc_trans.tet_id);
      voro_faces_ids.push_back(plane);
      voro_faces_sf_ids.push_back(hp.y == -1 ? hp.x : -1);
      voro_faces_neigh_id.push_back(
          hp.y == -1 ? -1 : (hp.x == cc_trans.voro_id ? hp.y : hp.x));
      voro_faces_num_adj.push_back(cc_trans.clip_trans(plane).h);
    }

    lf++;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// MainGuiWindow functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void MainGuiWindow::compute_partial_rpd_w_degeneration(
    std::vector<MedialSphere>& all_medial_spheres, const bool is_given_all_tets,
    bool is_debug) {
  instance_->rpd3d->calculate_partial(instance_->num_itr_global,
                                      instance_->num_sphere_added,
                                      is_given_all_tets);

  printf("--------------------Clean Degenerated Spheres--------------------\n");
  while (instance_->num_sphere_added = delete_degenerated_medial_spheres(
                                           all_medial_spheres, is_debug) > 0) {
    instance_->rpd3d->calculate_partial(instance_->num_itr_global,
                                        instance_->num_sphere_added,
                                        false /*is_given_all_tets*/);
  }
}
void MainGuiWindow::compute_rpd_init() {
  instance_->num_sphere_added = all_medial_spheres->size();
  compute_partial_rpd_w_degeneration(*(instance_->all_medial_spheres),
                                     true /*is_given_all_tets*/,
                                     false /*is_debug*/);

  // compute medial mesh
  instance_->run_mmesh_generator(*(instance_->sf_mesh), *(instance_->mmesh),
                                 *(instance_->all_medial_spheres),
                                 false /*is_trace_mstruc*/,
                                 false /*is_update_tan_pls*/);
}

void MainGuiWindow::show(bool is_compute_rpd) {
  // a few camera options
  polyscope::view::upDir = polyscope::UpDir::ZUp;

  // Initialize Polyscope
  polyscope::init();

  // Add the callback
  polyscope::state::userCallback = MainGuiWindow::callbacks;

  // show input tet_mesh
  show_tet_mesh(*(tet_mesh));

  // create RPD
  if (is_compute_rpd) compute_rpd_init();

  instance_->show_result_convex_cells(instance_->rpd3d->powercells,
                                      false /*is_slice_plane*/);
  instance_->show_medial_mesh(*(instance_->mmesh));
  instance_->show_medial_edges(*(instance_->mmesh));

  // Show the GUI
  polyscope::show();
}

void MainGuiWindow::auto_all() {
  timer::time_point start, end;
  start = timer::now();
  // Initialize Polyscope
  // polyscope::init();

  // Add the callback
  // polyscope::state::userCallback = MainGuiWindow::callbacks;

  // // show input tet_mesh
  // show_tet_mesh(*(tet_mesh));

  std::string name_no_ext = get_only_file_name(
      instance_->tet_mesh->tet_path_with_ext, false /*withExtension*/);
  // save_input_tet_houdini(*instance_->tet_mesh, name_no_ext);

  num_sphere_added = all_medial_spheres->size();
  instance_->rpd3d->calculate_partial(instance_->num_itr_global,
                                      instance_->num_sphere_added,
                                      true /*is_given_all_tets*/);

  auto_fix_extf();
  auto_fix_intf();
  auto_fix_extf();
  auto_fix_geo();
  auto_fix_intf();
  auto_fix_extf();
  auto_fix_topo();

  // save mat_dup_cnt
  // export_ma_dup_cnt(name_no_ext, *(instance_->mmesh));
  // save spheres， before thinning
  // save_spheres_file(*(instance_->all_medial_spheres), name_no_ext,
  //                   true /*is_save_type*/, true /*is_load_deleted*/,
  //                   true /*is_save_dup_cnt*/);

  // thinning
  prune(*(instance_->all_medial_spheres), *(instance_->mmesh),
        instance_->given_thinning_thres, true /*is_dup_cnt*/,
        false /*is_import_given*/, false /*is_sort_randomly*/,
        false /*is_debug*/);
  compute_Euler(*(instance_->mmesh));

  // save mat
  // export_ma(name_no_ext, *(instance_->mmesh));// not for computing Euler
  export_ma_clean(name_no_ext, *(instance_->mmesh));
  export_ma_ply(name_no_ext, *(instance_->mmesh));  // for computing Euler

  // // save rpd
  // save_convex_cells_houdini(*(instance_->params),
  //                           *(instance_->all_medial_spheres),
  //                           instance_->rpd3d->powercells, name_no_ext, false
  //                           /*is_slice_plane*/);
  end = timer::now();
  std::cout
      << "auto_all took "
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << " seconds in total" << std::endl;
}

void MainGuiWindow::auto_all_non_cad() {
  timer::time_point start, end;
  start = timer::now();
  instance_->num_sphere_added = all_medial_spheres->size();
  compute_partial_rpd_w_degeneration(*this->all_medial_spheres,
                                     true /*is_given_all_tets*/,
                                     false /*is_debug*/);

  auto_fix_geo();
  auto_fix_intf();
  auto_fix_extf();
  auto_fix_topo();

  std::string name_no_ext = get_only_file_name(
      instance_->tet_mesh->tet_path_with_ext, false /*withExtension*/);

  // save mat_dup_cnt
  // export_ma_dup_cnt(name_no_ext, *(instance_->mmesh));
  // save spheres， before thinning
  // save_spheres_file(*(instance_->all_medial_spheres), name_no_ext,
  //                   true /*is_save_type*/, true /*is_load_deleted*/,
  //                   true /*is_save_dup_cnt*/);

  // thinning
  instance_->given_thinning_thres = 0.1;
  prune(*(instance_->all_medial_spheres), *(instance_->mmesh),
        instance_->given_thinning_thres, true /*is_dup_cnt*/,
        false /*is_import_given*/, false /*is_sort_randomly*/,
        false /*is_debug*/);
  compute_Euler(*(instance_->mmesh));

  // save mat
  // export_ma(name_no_ext, *(instance_->mmesh));  // not for computing Euler
  export_ma_clean(name_no_ext, *(instance_->mmesh));  // for computing Euler
  export_ma_ply(name_no_ext, *(instance_->mmesh));    // for computing Euler

  end = timer::now();
  std::cout
      << "auto_all_non_cad took "
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << " seconds in total" << std::endl;
}

void MainGuiWindow::auto_fix_extf() {
  instance_->num_extf_itr = 0;
  int itr_max = this->is_non_cad ? EXTF_ITR_MAX : INT_MAX;
  timer::time_point start, end;
  start = timer::now();
  while (
      instance_->num_extf_itr < itr_max &&
      (instance_->num_sphere_added = instance_->run_extf_fix_itr(
           *(instance_->params), *(instance_->tet_mesh), *(instance_->sf_mesh),
           instance_->rpd3d->powercells, *(instance_->all_medial_spheres),
           *(instance_->mmesh), false /*is_debug*/)) > 0) {
    // instance_->rpd3d->calculate_partial(instance_->num_itr_global,
    //                                     instance_->num_sphere_added,
    //                                     false /*is_given_all_tets*/);
    instance_->compute_partial_rpd_w_degeneration(
        *(instance_->all_medial_spheres), false /*is_given_all_tets*/,
        false /*is_debug*/);
  }
  instance_->run_mmesh_generator(*(instance_->sf_mesh), *(instance_->mmesh),
                                 *(instance_->all_medial_spheres));
  end = timer::now();
  std::cout
      << "ExtfFix " << instance_->num_extf_itr << " took "
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << " seconds " << std::endl;
}

void MainGuiWindow::auto_fix_intf() {
  instance_->num_intf_itr = 0;
  int itr_max =
      this->is_non_cad ? instance_->num_intf_itr + INTF_ITR_MAX : INT_MAX;
  timer::time_point start, end;
  start = timer::now();
  while (
      instance_->num_intf_itr < itr_max &&
      (instance_->num_sphere_added = instance_->run_intf_fix_itr(
           *(instance_->sf_mesh), *(instance_->tet_mesh), *(instance_->mmesh),
           *(instance_->all_medial_spheres), false /*is_debug*/)) > 0) {
    // instance_->rpd3d->calculate_partial(instance_->num_itr_global,
    //                                     instance_->num_sphere_added,
    //                                     false /*is_given_all_tets*/);
    instance_->compute_partial_rpd_w_degeneration(
        *(instance_->all_medial_spheres), false /*is_given_all_tets*/,
        false /*is_debug*/);
    instance_->run_mmesh_generator(*(instance_->sf_mesh), *(instance_->mmesh),
                                   *(instance_->all_medial_spheres),
                                   false /*is_trace_mstruc*/);
  }  // while

  end = timer::now();
  std::cout
      << "IntfFix " << instance_->num_intf_itr << " took "
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << " seconds " << std::endl;
}

void MainGuiWindow::auto_fix_topo() {
  instance_->num_topo_itr = 0;
  int itr_max = this->is_non_cad ? instance_->num_topo_itr + 41 : TOPO_ITR_MAX;
  timer::time_point start, end;
  start = timer::now();
  // run topo fix
  while (instance_->num_topo_itr < itr_max &&
         (instance_->num_sphere_added = instance_->run_topo_fix_itr(
              *(instance_->sf_mesh), *(instance_->tet_mesh),
              *(instance_->all_medial_spheres), instance_->rpd3d->powercells,
              false /*is_debug*/)) > 0) {
    // instance_->rpd3d->calculate_partial(instance_->num_itr_global,
    //                                     instance_->num_sphere_added,
    //                                     false /*is_given_all_tets*/);
    instance_->compute_partial_rpd_w_degeneration(
        *(instance_->all_medial_spheres), false /*is_given_all_tets*/,
        false /*is_debug*/);
  }
  instance_->run_mmesh_generator(*(instance_->sf_mesh), *(instance_->mmesh),
                                 *(instance_->all_medial_spheres));
  end = timer::now();
  std::cout
      << "TopoFix " << instance_->num_topo_itr << " took "
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << " seconds " << std::endl;
}

void MainGuiWindow::auto_fix_geo() {
  instance_->num_geo_itr = 0;
  int itr_max = this->is_non_cad ? instance_->num_geo_itr + 5
                                 : instance_->num_geo_itr + GEO_ITR_MAX;
  timer::time_point start, end;
  start = timer::now();
  while (instance_->num_geo_itr < itr_max &&
         (instance_->num_sphere_added = instance_->run_fix_geo_itr(
              *(instance_->sf_mesh), *(instance_->tet_mesh),
              *(instance_->params), *(instance_->mmesh),
              *(instance_->all_medial_spheres), false /*is_debug*/)) > 0) {
    // instance_->rpd3d->calculate_partial(instance_->num_itr_global,
    //                                     instance_->num_sphere_added,
    //                                     false /*is_given_all_tets*/);
    instance_->compute_partial_rpd_w_degeneration(
        *(instance_->all_medial_spheres), false /*is_given_all_tets*/,
        false /*is_debug*/);
    instance_->run_mmesh_generator(*(instance_->sf_mesh), *(instance_->mmesh),
                                   *(instance_->all_medial_spheres));
  }
  end = timer::now();
  std::cout
      << "GeoFix " << instance_->num_geo_itr << " took "
      << std::chrono::duration_cast<std::chrono::seconds>(end - start).count()
      << " seconds " << std::endl;
}

// Your callback functions
void MainGuiWindow::callbacks() {
  ImGui::PushItemWidth(100);

  // Surface
  if (ImGui::Button("show surface mesh")) {
    instance_->show_surface_mesh(*(instance_->sf_mesh),
                                 -1 /*given_sf_face_id*/);
  }
  // Medial spheres
  if (ImGui::SmallButton("show all mat centers")) {
    instance_->show_all_medial_spheres(*(instance_->all_medial_spheres));
  }

  ImGui::PopItemWidth();
}

void MainGuiWindow::show_result_convex_cells(
    std::vector<ConvexCellHost>& powercells, bool is_slice_plane) {
  std::vector<cfloat3> voro_points;
  std::vector<float> voro_points_num_adjs;
  std::vector<std::array<int, 4>> voro_tets;  // tet
  std::vector<adouble3> voro_colors;          // random color per site
  std::vector<float> voro_tets_euler;
  std::vector<int> voro_tets_sites, voro_tets_cell_ids;
  std::map<int, adouble3> color_map_per_site;

  int num_old_tet = 0;
  for (auto& cell_trans : powercells) {
    assert(cell_trans.id != -1);
    // // filter spheres by slicing plane
    // Vector3 last_bary = compute_cell_barycenter(cell_trans);
    // if (is_slice_plane && is_slice_by_plane(last_bary, *(this->params)))
    //   continue;

    // if (cell_trans.voro_id != 5 && cell_trans.voro_id != 19 &&
    //     cell_trans.voro_id != 29)
    //   continue;
    num_old_tet = voro_tets.size();
    get_convex_cell_tet(cell_trans, voro_points, voro_tets, voro_tets_sites,
                        voro_tets_cell_ids, voro_tets_euler);
    // assign random color per site
    if (color_map_per_site.find(cell_trans.voro_id) == color_map_per_site.end())
      color_map_per_site[cell_trans.voro_id] = {{polyscope::randomUnit(),
                                                 polyscope::randomUnit(),
                                                 polyscope::randomUnit()}};
    FOR(i, voro_tets.size() - num_old_tet) {
      voro_colors.push_back(color_map_per_site.at(cell_trans.voro_id));
    }
  }

  // Register the volume mesh with Polyscope
  auto my_mesh = polyscope::registerTetMesh("my RPD", voro_points, voro_tets);
  // Add a scalar quantity
  my_mesh->addCellScalarQuantity("site", voro_tets_sites);
  my_mesh->addCellColorQuantity("site_color", voro_colors)->setEnabled(true);
  my_mesh->addCellScalarQuantity("euler", voro_tets_euler);
  my_mesh->addCellScalarQuantity("cell_id", voro_tets_cell_ids);
}

void MainGuiWindow::clear_site_boundary() {
  for (int i = 0; i < 100000; i++) {
    std::string site_name = "site" + std::to_string(i);
    polyscope::removeSurfaceMesh(site_name, false);
  }
}

void MainGuiWindow::show_one_site_boundary(const int given_site_id) {
  auto& all_cells = this->rpd3d->powercells;
  int cnt = 0;
  // assign site a color
  adouble3 ran_color = {{polyscope::randomUnit(), polyscope::randomUnit(),
                         polyscope::randomUnit()}};

  // variables
  std::vector<cfloat3> voro_points;
  std::vector<float> voro_points_num_adjs;
  std::vector<aint3> voro_faces;  // face
  std::vector<float> voro_faces_cell_eulers;
  std::vector<int> voro_faces_sites, voro_faces_cell_ids, voro_faces_tet_ids,
      voro_faces_lids, voro_faces_num_adjs, voro_faces_sf_ids,
      voro_faces_neigh_id;

  // loop all convex cells
  for (auto& one_cell_trans : all_cells) {
    if (one_cell_trans.voro_id != given_site_id) continue;
    get_convex_cell_faces(
        one_cell_trans, voro_points, voro_points_num_adjs, voro_faces,
        voro_faces_sites, voro_faces_cell_ids, voro_faces_tet_ids,
        voro_faces_cell_eulers, voro_faces_lids, voro_faces_sf_ids,
        voro_faces_neigh_id, voro_faces_num_adjs, this->sf_mesh->facets.nb(),
        true /*is_boundary_only*/);
  }

  // Register the surface mesh with Polyscope
  std::string site_name = "site" + std::to_string(given_site_id);
  // std::string site_name = "site";
  auto my_mesh =
      polyscope::registerSurfaceMesh(site_name, voro_points, voro_faces);
  my_mesh->addVertexScalarQuantity("local num_adjs", voro_points_num_adjs);
  my_mesh->addFaceScalarQuantity("site", voro_faces_sites);
  my_mesh->addFaceScalarQuantity("cell euler", voro_faces_cell_eulers);
  my_mesh->addFaceScalarQuantity("local fid", voro_faces_lids);
  my_mesh->addFaceScalarQuantity("sf_mesh fid", voro_faces_sf_ids);
  my_mesh->addFaceScalarQuantity("neigh site id", voro_faces_neigh_id);
  my_mesh->addFaceScalarQuantity("local num_adjs", voro_faces_num_adjs);
  my_mesh->addFaceScalarQuantity("cell id", voro_faces_cell_ids);
  my_mesh->addFaceScalarQuantity("tet id", voro_faces_tet_ids);
  std::vector<adouble3> voro_colors;
  for (int i = 0; i < voro_faces.size(); i++) voro_colors.push_back(ran_color);
  my_mesh->addFaceColorQuantity("color", voro_colors)->setEnabled(true);
}

void MainGuiWindow::clear_cells() {
  std::string cell_name = "cell";
  for (uint i = 0; i < 3000; i++) {
    for (uint j = 0; j < this->site_cnt; j++) {
      std::string new_cell_name =
          cell_name + std::to_string(j) + "_" + std::to_string(i);
      polyscope::removeSurfaceMesh(new_cell_name, false);
    }
  }
  this->site_cnt = 0;
}

void MainGuiWindow::show_one_site(const int given_site_id) {
  auto& all_cells = this->rpd3d->powercells;
  int cnt = 0;
  // assign site a color
  adouble3 ran_color = {{polyscope::randomUnit(), polyscope::randomUnit(),
                         polyscope::randomUnit()}};
  for (const auto& one_cell_trans : all_cells) {
    if (one_cell_trans.voro_id != given_site_id) continue;
    // if (one_cell_trans.id != 3303 && one_cell_trans.id != 2031) continue;
    // printf("showing site %d cell %d \n", given_site_id, one_cell_trans.id);
    //
    show_one_cell(one_cell_trans.id, ran_color, cnt++);
  }
}

void MainGuiWindow::show_one_cell(const int given_cell_id, adouble3 ran_color,
                                  int suffix) {
  auto& all_cells = this->rpd3d->powercells;
  if (given_cell_id < 0 || given_cell_id >= all_cells.size()) return;
  auto& cell_trans = all_cells[given_cell_id];
  bool is_color_given = false;
  if (ran_color[0] != -1 && ran_color[1] != -1 && ran_color[2] != -1)
    is_color_given = true;

  std::vector<cfloat3> voro_points;
  std::vector<float> voro_points_num_adjs;
  std::vector<aint3> voro_faces;  // face
  std::vector<float> voro_faces_cell_eulers;
  std::vector<int> voro_faces_sites, voro_faces_cell_ids, voro_faces_tet_ids,
      voro_faces_lids, voro_faces_num_adjs, voro_faces_sf_ids,
      voro_faces_neigh_id;

  get_convex_cell_faces(cell_trans, voro_points, voro_points_num_adjs,
                        voro_faces, voro_faces_sites, voro_faces_cell_ids,
                        voro_faces_tet_ids, voro_faces_cell_eulers,
                        voro_faces_lids, voro_faces_sf_ids, voro_faces_neigh_id,
                        voro_faces_num_adjs, this->sf_mesh->facets.nb(),
                        false /*is_boundary_only*/);

  // printf("cell %d has points: \n", given_cell_id);
  // for (const auto& p : voro_points) {
  //   printf("v %f %f %f \n", p.x, p.y, p.z);
  // }
  // for (const auto& f : voro_faces) {
  //   printf("f %d %d %d\n", f[0], f[1], f[2]);
  // }

  // Register the surface mesh with Polyscope
  std::string cell_name = "cell";
  // std::string cell_name = "cell" + std::to_string(given_cell_id);
  if (suffix != -1) {
    cell_name += std::to_string(this->site_cnt) + "_" + std::to_string(suffix);
    this->site_cnt++;
  }
  auto my_mesh =
      polyscope::registerSurfaceMesh(cell_name, voro_points, voro_faces);
  my_mesh->addVertexScalarQuantity("local num_adjs", voro_points_num_adjs);
  my_mesh->addFaceScalarQuantity("site", voro_faces_sites);
  my_mesh->addFaceScalarQuantity("cell euler", voro_faces_cell_eulers);
  my_mesh->addFaceScalarQuantity("local fid", voro_faces_lids);
  my_mesh->addFaceScalarQuantity("sf_mesh fid", voro_faces_sf_ids);
  my_mesh->addFaceScalarQuantity("neigh site id", voro_faces_neigh_id);
  my_mesh->addFaceScalarQuantity("local num_adjs", voro_faces_num_adjs);
  my_mesh->addFaceScalarQuantity("cell id", voro_faces_cell_ids);
  my_mesh->addFaceScalarQuantity("tet id", voro_faces_tet_ids);
  if (is_color_given) {
    std::vector<adouble3> voro_colors;
    for (int i = 0; i < voro_faces.size(); i++)
      voro_colors.push_back(ran_color);
    my_mesh->addFaceColorQuantity("color", voro_colors)->setEnabled(true);
  }
}

void MainGuiWindow::show_tet_mesh(const TetMesh& tet_mesh) {
  std::vector<std::array<float, 3>> tet_vertices_new;
  std::vector<std::array<int, 4>> tet_indices_new;
  std::vector<int> tet_ids;
  for (int i = 0; i < tet_mesh.tet_vertices.size() / 3; i++) {
    tet_vertices_new.push_back(
        {{tet_mesh.tet_vertices[i * 3], tet_mesh.tet_vertices[i * 3 + 1],
          tet_mesh.tet_vertices[i * 3 + 2]}});
  }
  for (int i = 0; i < tet_mesh.tet_indices.size() / 4; i++) {
    tet_indices_new.push_back(
        {{tet_mesh.tet_indices[i * 4], tet_mesh.tet_indices[i * 4 + 1],
          tet_mesh.tet_indices[i * 4 + 2], tet_mesh.tet_indices[i * 4 + 3]}});
    tet_ids.push_back(i);
  }
  auto my_tet =
      polyscope::registerTetMesh("my tet", tet_vertices_new, tet_indices_new);
  my_tet->addCellScalarQuantity("tet_id", tet_ids);
  my_tet->setEnabled(false);
}

void MainGuiWindow::show_one_tet(const TetMesh& tet_mesh,
                                 const int& given_tet_id) {
  int nb_tet = tet_mesh.tet_indices.size() / 4;
  if (given_tet_id < 0 || given_tet_id >= nb_tet) {
    polyscope::warning("tet range is [0, " + std::to_string(nb_tet - 1) + "]");
    return;
  }
  std::vector<std::array<float, 3>> tet_vertices_new;
  std::vector<std::array<int, 4>> tet_indices_new;
  std::vector<int> tet_ids;
  tet_indices_new.push_back({{0, 1, 2, 3}});
  tet_ids.push_back(given_tet_id);

  for (int i = 0; i < 4; i++) {
    std::array<float, 3> v;
    for (int j = 0; j < 3; j++) {
      int idx = tet_mesh.tet_indices[given_tet_id * 4 + i];
      v[j] = tet_mesh.tet_vertices[idx * 3 + j];
    }
    tet_vertices_new.push_back(v);
  }

  auto one_tet =
      polyscope::registerTetMesh("one tet", tet_vertices_new, tet_indices_new);
  one_tet->addCellScalarQuantity("tet_id", tet_ids);
}

void MainGuiWindow::show_surface_mesh(const GEO::Mesh& sf_mesh,
                                      const int& given_sf_face_id) {
  std::vector<adouble3> sf_vertices;
  std::vector<aint3> sf_faces;  // vs size might be > 3
  std::vector<bool> is_selected(sf_mesh.facets.nb(), false);
  std::vector<int> sf_vs_tet_id_attr;
  const GEO::Attribute<int> tet_vid_attr(sf_mesh.vertices.attributes(),
                                         "tet_vid");

  std::vector<Vector3> corner_vertices;
  std::set<aint2> s_edges, cc_edges;
  GEO::Attribute<int> attr_corners(sf_mesh.vertices.attributes(), "corner");
  GEO::Attribute<int> attr_se(sf_mesh.edges.attributes(), "se");
  GEO::Attribute<int> attr_cce(sf_mesh.edges.attributes(), "cce");

  for (uint v = 0; v < sf_mesh.vertices.nb(); v++) {
    const Vector3& p = sf_mesh.vertices.point(v);
    sf_vertices.push_back({{p[0], p[1], p[2]}});
    sf_vs_tet_id_attr.push_back(tet_vid_attr[v]);
    if (attr_corners[v] == 1) corner_vertices.push_back(p);
  }

  for (uint f = 0; f < sf_mesh.facets.nb(); f++) {
    if (given_sf_face_id == f) {
      is_selected[f] = true;
    }
    int f_nb_v = sf_mesh.facets.nb_vertices(f);
    assert(f_nb_v == 3);
    aint3 one_f;
    for (uint lv = 0; lv < f_nb_v; lv++) {
      int v = sf_mesh.facets.vertex(f, lv);
      one_f[lv] = v;
    }
    sf_faces.push_back(one_f);
  }

  for (uint e = 0; e < sf_mesh.edges.nb(); e++) {
    aint2 edge;
    for (int lv = 0; lv < 2; ++lv) {
      edge[lv] = sf_mesh.edges.vertex(e, lv);
    }
    if (attr_se[e] == 1) {
      s_edges.insert(edge);
    } else if (attr_cce[e] == 1) {
      cc_edges.insert(edge);
    }
  }

  auto poly_sf_mesh =
      polyscope::registerSurfaceMesh("Surface mesh", sf_vertices, sf_faces);
  poly_sf_mesh->addFaceScalarQuantity("is_selected", is_selected)
      ->setEnabled(true);
  poly_sf_mesh->addVertexScalarQuantity("tet_vid", sf_vs_tet_id_attr);
  auto poly_corners = polyscope::registerPointCloud("Corners", corner_vertices);
  poly_corners->setPointRadius(MainGuiWindow::point_radius_rel);
  convert_to_show_edges(sf_vertices, s_edges, "Sharp Edges");
  convert_to_show_edges(sf_vertices, cc_edges, "Concave Edges");
}

void MainGuiWindow::convert_to_show_edges(const std::vector<adouble3>& old_pos,
                                          const std::set<aint2>& old_edges,
                                          std::string curve_name) {
  if (old_pos.empty() || old_edges.empty()) {
    // todo: remove existing edge
    return;
  }

  int new_idx = 0;
  std::map<int, int> old2new;
  std::vector<adouble3> new_pos;
  for (const auto& one_old_edge : old_edges) {
    for (const auto& v : one_old_edge) {
      if (old2new.find(v) != old2new.end()) {
        continue;
      }
      old2new[v] = new_idx;
      new_pos.push_back(old_pos[v]);
      new_idx++;
    }
  }

  std::set<aint2> new_edges;
  for (const auto& one_old_edge : old_edges) {
    aint2 one_new_edge = {{old2new[one_old_edge[0]], old2new[one_old_edge[1]]}};
    std::sort(one_new_edge.begin(), one_new_edge.end());
    new_edges.insert(one_new_edge);
  }

  auto poly_edges =
      polyscope::registerCurveNetwork(curve_name, new_pos, new_edges);
  poly_edges->setRadius(MainGuiWindow::edge_radius_rel, true /*isRelative*/);
}

void MainGuiWindow::show_pin_points(
    const std::vector<MedialSphere>& all_medial_spheres,
    const int given_sphere_id) {
  if (given_sphere_id == -1) return;
  std::vector<Vector3> all_pins;
  std::vector<Vector3> all_pin_normals;
  for (uint i = 0; i < all_medial_spheres.size(); i++) {
    if (given_sphere_id != -1 && i != given_sphere_id) continue;
    const auto& msphere = all_medial_spheres.at(i);
    all_pins.push_back(msphere.ss.p);
    all_pin_normals.push_back(msphere.ss.p_normal);
  }

  auto all_pin = polyscope::registerPointCloud("All pins", all_pins);
  all_pin->setEnabled(false);
  all_pin->addVectorQuantity("normal", all_pin_normals)->setEnabled(false);
}

void MainGuiWindow::show_all_medial_spheres(
    const std::vector<MedialSphere>& all_medial_spheres) {
  std::vector<Vector3> all_spheres;
  std::vector<int> sphere_ids, is_on_ce_pin, is_radius_dilated, sphere_types,
      sphere_itr_cnt;
  for (uint i = 0; i < all_medial_spheres.size(); i++) {
    const auto& msphere = all_medial_spheres.at(i);
    if (msphere.is_deleted) continue;
    all_spheres.push_back(msphere.center);
    sphere_ids.push_back(msphere.id);
    is_on_ce_pin.push_back(msphere.is_on_ce_pin() ? 1 : 0);
    is_radius_dilated.push_back(msphere.is_radius_dilated ? 1 : 0);
    sphere_types.push_back(msphere.type);
    sphere_itr_cnt.push_back(msphere.itr_cnt);
  }
  auto all_pc = polyscope::registerPointCloud("all spheres", all_spheres);
  all_pc->setPointRadius(MainGuiWindow::point_radius_rel)->setEnabled(true);
  all_pc->addScalarQuantity("id", sphere_ids);
  all_pc->addScalarQuantity("is_on_ce_pin", is_on_ce_pin)->setEnabled(false);
  all_pc->addScalarQuantity("is_radius_dilated", is_radius_dilated);
  all_pc->addScalarQuantity("type", sphere_types);
  all_pc->addScalarQuantity("itr_cnt", sphere_itr_cnt)->setEnabled(true);
}

void MainGuiWindow::show_one_sphere(
    const std::vector<MedialSphere>& all_medial_spheres,
    const int given_sphere_id, const bool is_show_radius,
    const bool is_clear_all) {
  if (given_sphere_id < 0 || given_sphere_id >= all_medial_spheres.size()) {
    polyscope::warning("msphere tag range is [0, " +
                       std::to_string(all_medial_spheres.size() - 1) + "]");
    return;
  }
  std::vector<Vector3> positions;
  std::vector<double> radii;
  const auto& msphere = all_medial_spheres.at(given_sphere_id);
  positions.push_back(msphere.center);
  if (msphere.is_on_se())
    radii.push_back(point_radius_rel);
  else
    radii.push_back(msphere.radius);
  msphere.print_info();

  std::string name_sphere = "Sphere ";
  // if is_clear_all, then shows nothing
  if (is_clear_all) {
    for (uint i = 1; i < 20; i++) {
      std::string new_sphere =
          name_sphere + "(with radius)" + std::to_string(i);
      polyscope::removePointCloud(new_sphere, false);
      new_sphere = name_sphere + "(no radius)" + std::to_string(i);
      polyscope::removePointCloud(new_sphere, false);
    }
  }

  if (is_show_radius)
    name_sphere += "(with radius)";
  else
    name_sphere += "(no radius)";
  for (uint i = 1; i < 20; i++) {
    std::string new_sphere = name_sphere + std::to_string(i);
    if (polyscope::hasPointCloud(new_sphere) == false) {
      auto single_msphere =
          polyscope::registerPointCloud(new_sphere, positions);
      if (is_show_radius) {
        single_msphere->addScalarQuantity("radius", radii);
        single_msphere->setPointRadiusQuantity("radius", false);
      }
      break;
    }
  }

  // Tangent planes
  std::vector<Vector3> tan_pl_normals, tan_pins;
  for (const auto& tan_pl : msphere.tan_planes) {
    tan_pl_normals.push_back(tan_pl.normal);
    tan_pins.push_back(tan_pl.tan_point);
  }
  // Tangent concave lines
  for (const auto& tan_cc_line : msphere.tan_cc_lines) {
    tan_pl_normals.push_back(tan_cc_line.normal);
    tan_pins.push_back(tan_cc_line.tan_point);
  }
  auto mat_p_tan_elements =
      polyscope::registerPointCloud("Tangent Elements", tan_pins);
  mat_p_tan_elements->addVectorQuantity("normal", tan_pl_normals)
      ->setEnabled(true);
}

void MainGuiWindow::show_medial_mesh(const MedialMesh& mmesh,
                                     int given_mat_face_id) {
  if (mmesh.vertices == nullptr) return;
  const auto& mspheres = *(mmesh.vertices);
  const auto& mfaces = mmesh.faces;

  std::vector<Vector3> mat_pos(mspheres.size());
  std::vector<double> mat_radius(mspheres.size());
  std::vector<aint3> mat_faces(mfaces.size(), {{0, 0, 0}});
  std::vector<double> mat_faces_importance(mfaces.size(), 2.0);
  std::vector<int> mat_faces_cnt(mfaces.size(), 0);

  // set true when given_mat_face_id != -1
  std::vector<bool> is_given_mat_face(mmesh.faces.size(), false);
  // std::vector<int> mat_faces_mstruc_ids(mfaces.size(), -1);
  std::vector<Vector3> dual_segment_vs;
  std::vector<aint2> dual_segment_edge;
  // for simple triangles
  std::vector<Vector3> st_points;
  std::vector<aint3> st_faces;
  std::vector<Vector3> st_centroids;
  std::vector<Vector3> st_normals;
  std::vector<Vector3> st_nearest_p;

  // store mat vertices
  for (uint i = 0; i < mspheres.size(); i++) {
    mat_pos[i] = mspheres[i].center;
    mat_radius[i] = mspheres[i].radius;
  }

  // store mat faces
  // need to enbale culling for drawing ma faces twice
  for (uint f = 0; f < mfaces.size(); f++) {
    if (mfaces[f].is_deleted) continue;
    // if (mfaces[f].dup_cnt <= 1) continue;
    mat_faces_cnt[f] = mfaces[f].dup_cnt;
    mat_faces_importance[f] = mfaces[f].importance;
    // mat_faces_mstruc_ids[f] = mfaces[f].mstruc_id;
    // draw facets
    uint lv = 0;
    for (const auto& v : mfaces[f].vertices_) {
      mat_faces[f][lv++] = v;
    }
  }

  // if this mat face is what we are looking for
  if (given_mat_face_id > 0 && given_mat_face_id < mmesh.faces.size()) {
    int f = given_mat_face_id;
    const auto& matf = mmesh.faces[f];
    is_given_mat_face[f] = true;
    for (int i = 0; i < matf.dual_edge_endpoints.size(); i++) {
      dual_segment_vs.push_back(matf.dual_edge_endpoints[i][0].first);
      dual_segment_vs.push_back(matf.dual_edge_endpoints[i][1].first);
      dual_segment_edge.push_back({{i * 2, i * 2 + 1}});
    }
    // show simple triangles
    for (int i = 0; i < 2; i++) {
      st_centroids.push_back(matf.st[i].centroid);
      st_normals.push_back(matf.st[i].normal);
      st_nearest_p.push_back(matf.st[i].nearest_point);
      for (int j = 0; j < 3; j++) st_points.push_back(matf.st[i].v[j]);
      st_faces.push_back({{i * 3, i * 3 + 1, i * 3 + 2}});
    }

    matf.print_medial_face();
    for (const auto tid : matf.tets_) {
      mmesh.tets.at(tid).print_medial_tet();
    }
    for (int i = 0; i < 2; i++) {
      printf("mface %d has nearest surface fid %d, dist_to_sf %f\n", matf.fid,
             matf.st[i].nearest_sf_fid, matf.st[i].dist_to_sf);
    }
  }

  // Register medial mesh
  auto medial_mesh =
      polyscope::registerSurfaceMesh("my medial mesh", mat_pos, mat_faces);
  medial_mesh->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
  medial_mesh->addVertexScalarQuantity("radius", mat_radius,
                                       polyscope::DataType::MAGNITUDE);
  medial_mesh->addFaceScalarQuantity("importance", mat_faces_importance)
      ->setEnabled(false);
  // medial_mesh->addFaceScalarQuantity("dup_cnt", mat_faces_cnt)
  //     ->setEnabled(false);
  // medial_mesh->addFaceScalarQuantity("mstruc_id", mat_faces_mstruc_ids)
  //     ->setEnabled(true);

  if (given_mat_face_id != -1) {
    medial_mesh->addFaceScalarQuantity("is_given_mat_face", is_given_mat_face)
        ->setEnabled(true);
    polyscope::registerCurveNetwork("dual_segment", dual_segment_vs,
                                    dual_segment_edge);
    // polyscope::registerSurfaceMesh("mface_st", st_points, st_faces);
    // auto tmp = polyscope::registerPointCloud("st centroids", st_centroids);
    // tmp->addVectorQuantity("st normals", st_normals)->setEnabled(true);
    // polyscope::registerPointCloud("mface_st nearest_p", st_nearest_p);
  }
}

void MainGuiWindow::show_medial_edges(const MedialMesh& mmesh) {
  const auto& mspheres = *(mmesh.vertices);
  const auto& medges = mmesh.edges;
  std::vector<Vector3> mat_pos(mspheres.size());
  std::vector<aint2> mat_edges;
  std::vector<int> mat_edges_extf;
  // store mat vertices
  for (uint i = 0; i < mspheres.size(); i++) {
    mat_pos[i] = mspheres[i].center;
  }

  // store mat edges
  for (uint e = 0; e < medges.size(); e++) {
    if (medges[e].is_deleted) continue;
    mat_edges.push_back(medges[e].vertices_);
    mat_edges_extf.push_back(medges[e].is_extf);
  }

  auto me = polyscope::registerCurveNetwork("medial edges", mat_pos, mat_edges);
  me->setRadius(MainGuiWindow::edge_radius_rel, true /*isRelative*/);
  me->addEdgeScalarQuantity("extf", mat_edges_extf)->setEnabled(true);
}

void MainGuiWindow::show_medial_edge_info(MedialMesh& mmesh,
                                          int given_mat_edge_v1,
                                          int given_mat_edge_v2) {
  if (given_mat_edge_v1 < 0 || given_mat_edge_v1 >= mmesh.vertices->size()) {
    printf("given_mat_edge_v1 %d out of bound [0, %zu)\n", given_mat_edge_v1,
           mmesh.edges.size());
    return;
  }
  if (given_mat_edge_v2 < 0 || given_mat_edge_v2 >= mmesh.vertices->size()) {
    printf("given_mat_edge_v2 %d out of bound [0, %zu)\n", given_mat_edge_v2,
           mmesh.edges.size());
    return;
  }
  int eid = -1;
  mmesh.Edge(given_mat_edge_v1, given_mat_edge_v2, eid);
  if (eid == -1) {
    printf("cannot find mmesh edge for (%d,%d)\n", given_mat_edge_v1,
           given_mat_edge_v2);
    return;
  }

  const auto& medge = mmesh.edges.at(eid);
  printf("medge %d (%d,%d), is_intf: %d\n", medge.is_intf);

  printf("medge %d (%d,%d), common_tan_pls: ----------------------------\n",
         eid, medge.vertices_[0], medge.vertices_[1]);
  for (const auto& common : medge.common_tan_pls) {
    common.print_info();
  }

  printf(
      "medge %d (%d,%d), v1_non_common_tan_pls: "
      "----------------------------\n",
      eid, medge.vertices_[0], medge.vertices_[1]);
  for (const auto& v1_non_common : medge.v1_non_common_tan_pls) {
    v1_non_common.print_info();
  }

  printf(
      "medge %d (%d,%d), v2_non_common_tan_pls: "
      "----------------------------\n",
      eid, medge.vertices_[0], medge.vertices_[1]);
  for (const auto& v_non_common : medge.v2_non_common_tan_pls) {
    v_non_common.print_info();
  }
}

void MainGuiWindow::show_fix_geo_samples(
    const std::vector<double>& geo_samples,
    const std::vector<float>& samples_dist2mat,
    const std::vector<aint3>& samples_clostprim, const int given_sample_id) {
  std::vector<Vector3> all_samples;
  int sample_size = geo_samples.size() / 3;
  for (uint i = 0; i < sample_size; i++) {
    if (given_sample_id != -1 && i != given_sample_id) continue;
    if (i == given_sample_id)
      printf("sample %d has dist2mat %f prim (%d,%d,%d)\n", i,
             samples_dist2mat.at(i), samples_clostprim.at(i)[0],
             samples_clostprim.at(i)[1], samples_clostprim.at(i)[2]);
    all_samples.push_back(Vector3(geo_samples.at(i * 3),
                                  geo_samples.at(i * 3 + 1),
                                  geo_samples.at(i * 3 + 2)));
  }

  auto all_pin = polyscope::registerPointCloud("Geo samples", all_samples);
}