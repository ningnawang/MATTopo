
#include "io_api.h"

#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <igl/writePLY.h>

#include <Eigen/Core>
#include <bitset>

#include "../extern/mshloader/MshLoader.h"
#include "io_utils.hpp"
#include "voronoi_common.h"

void load_spheres_from_file(const char* filename,
                            std::vector<MedialSphere>& all_medial_spheres,
                            bool is_load_type, bool is_load_deleted,
                            bool is_save_dup_cnt) {
  std::ifstream file(filename);
  int dim, n_site;
  int type, is_deleted, dup_cnt;
  file >> dim >> n_site;
  assert(dim == 3 || dim == 4);
  std::cout << "n_site: " << n_site << ", dim: " << dim << std::endl;

  all_medial_spheres.clear();
  for (int i = 0; i < n_site; i++) {
    MedialSphere msphere(all_medial_spheres.size(), Vector3(0, 0, 0),
                         Vector3(0, 0, 0));
    file >> msphere.center[0] >> msphere.center[1] >> msphere.center[2];
    if (dim == 4)
      file >> msphere.radius;
    else
      msphere.radius = 1;
    if (is_load_type) {
      file >> type;
      if (type != SphereType::T_UNK)
        msphere.type = SphereType(type);  // else as T2 sphere
    }
    if (is_load_deleted) {
      file >> is_deleted;
      msphere.is_deleted = false;
      if (is_deleted) msphere.is_deleted = true;
    }
    if (is_save_dup_cnt) {
      file >> dup_cnt;
      msphere.dup_cnt = dup_cnt;
    }
    all_medial_spheres.push_back(msphere);
    // printf("site %d has sphere: (%lf %lf %lf %lf), type: %d\n", i,
    //        msphere.center[0], msphere.center[1], msphere.center[2],
    //        msphere.radius, msphere.type);
  }
  file.close();
}

void save_spheres_file(const std::vector<MedialSphere>& all_medial_spheres,
                       const std::string filename, bool is_save_type,
                       bool is_load_deleted, bool is_save_dup_cnt) {
  std::string sphere_path =
      "../out/sph/sph_" + filename + "_" + get_timestamp() + ".sph";
  int n_site = all_medial_spheres.size();
  std::fstream file;
  file.open(sphere_path, std::ios_base::out);
  file << 4 << " " << n_site << std::endl;
  for (int i = 0; i < n_site; i++) {
    const auto& msphere = all_medial_spheres.at(i);
    file << std::setiosflags(std::ios::fixed) << std::setprecision(15)
         << msphere.center[0] << " " << msphere.center[1] << " "
         << msphere.center[2] << " " << msphere.radius;
    if (is_save_type) file << " " << msphere.type;
    if (is_load_deleted) file << " " << msphere.is_deleted;
    if (is_save_dup_cnt) file << " " << msphere.dup_cnt;
    file << std::endl;
  }
  file.close();
  printf("saved .sph file %s\n", sphere_path.c_str());
}

// helper function for export_ma_clean() and export_ma_ply()
void get_mat_clean(const MedialMesh& mat, std::vector<Vector4>& vertices,
                   std::vector<aint2>& edges,
                   std::vector<std::array<int, 3>>& faces) {
  std::map<int, int> map_vertices;  // mat vertex tag to new
  vertices.clear();
  edges.clear();
  faces.clear();

  auto get_vertex_mapped_id = [&](const int vid) {
    if (map_vertices.find(vid) == map_vertices.end()) {
      // add a new vertex
      map_vertices[vid] = map_vertices.size();
    }
    return map_vertices.at(vid);
  };

  // faces
  for (int f = 0; f < mat.faces.size(); f++) {
    const auto& face = mat.faces[f];
    if (face.is_deleted) continue;
    std::array<int, 3> one_f;
    for (uint j = 0; j < 3; j++) {
      int vid = get_vertex_mapped_id(face.vertices_[j]);
      one_f[j] = vid;
    }
    faces.push_back(one_f);
  }
  printf("faces: %d \n", faces.size());

  // edges
  for (int e = 0; e < mat.edges.size(); e++) {
    const auto& edge = mat.edges[e];
    if (edge.is_deleted) continue;
    // if (edge.faces_.empty()) continue;
    int vid1 = get_vertex_mapped_id(edge.vertices_[0]);
    int vid2 = get_vertex_mapped_id(edge.vertices_[1]);
    edges.push_back({{vid1, vid2}});
  }
  printf("edges: %d \n", edges.size());

  // vertices
  // save from map_vertices, to avoid (0,0,0,0) in .ma file
  vertices.resize(map_vertices.size());
  for (const auto& v_pair : map_vertices) {
    int old_vid = v_pair.first;
    int new_vid = v_pair.second;
    const auto& mat_v = mat.vertices->at(old_vid);
    vertices[new_vid] = Vector4(mat_v.center[0], mat_v.center[1],
                                mat_v.center[2], mat_v.radius);
  }
  // for (int v = 0; v < mat.vertices->size(); v++) {
  //   const auto& mat_v = mat.vertices->at(v);
  //   if (mat_v.is_deleted) continue;
  //   if (mat_v.edges_.empty() && mat_v.faces_.empty()) continue;  //
  //   isolated？ int vid = get_vertex_mapped_id(mat_v.id); if (vid ==
  //   map_vertices.size()) {
  //     // vertices.push_back(Vector4(vertex.center[0], vertex.center[1],
  //     //                            vertex.center[2], vertex.radius));
  //     printf("mat sphere: %d not connect to any edge/face, wrong?\n", v);
  //     assert(false);
  //   } else {
  //     vertices[vid] = Vector4(mat_v.center[0], mat_v.center[1],
  //     mat_v.center[2],
  //                             mat_v.radius);
  //   }
  // }
  printf("vertcies: %d \n", vertices.size());
}

// remove deleted edges/faces
// should save the same result as export_ma_ply but with different format
//
// this export is matching blender addon
// https://github.com/songshibo/blender-mat-addon
/* # .ma format
 * numVertices numEdges numFaces
 * v x y z r
 * e v1 v2
 * f v1 v2 v3
 */
void export_ma_clean(const std::string& maname, const MedialMesh& mat) {
  std::string ma_name_full =
      "../out/mat/mat_" + maname + "_" + get_timestamp() + ".ma";
  printf("start saving mat .m file: %s \n", ma_name_full);

  std::vector<Vector4> vertices;
  std::vector<aint2> edges;
  std::vector<std::array<int, 3>> faces;
  get_mat_clean(mat, vertices, edges, faces);

  // export
  export_ma_given(maname, vertices, edges, faces);
}

// remove deleted edges/faces
void export_ma_ply(const std::string& maname, const MedialMesh& mat) {
  std::string ma_name_full =
      "../out/mat/mat_" + maname + "_" + get_timestamp() + ".ply";
  // std::string ma_name_full = "../out/mat/mat_" + maname + ".ply";
  printf("start saving mat .ply file: %s \n", ma_name_full);

  std::vector<Vector4> vertices;
  std::vector<aint2> edges;
  std::vector<std::array<int, 3>> faces;
  get_mat_clean(mat, vertices, edges, faces);

  Eigen::MatrixXd V(vertices.size(), 3);
  Eigen::MatrixXi E(edges.size(), 2);
  Eigen::MatrixXi F(faces.size(), 3);
  for (int v = 0; v < vertices.size(); v++) {
    V(v, 0) = vertices[v][0];
    V(v, 1) = vertices[v][1];
    V(v, 2) = vertices[v][2];
  }
  for (int e = 0; e < edges.size(); e++) {
    E(e, 0) = edges[e][0];
    E(e, 1) = edges[e][1];
  }
  for (int f = 0; f < faces.size(); f++) {
    F(f, 0) = faces[f][0];
    F(f, 1) = faces[f][1];
    F(f, 2) = faces[f][2];
  }
  igl::writePLY(ma_name_full, V, F, E);
}

/*
 * numVertices numEdges numFaces numTets
 * v x y z r flag_type flag_delete dup_cnt
 * e v1 v2 flag_type flag_delete dup_cnt
 * f v1 v2 v3 importance flag_delete dup_cnt
 * t v1 v2 v3 v4 flag_delete dup_cnt
 */
void export_ma_dup_cnt(const std::string& maname, MedialMesh& mat) {
  std::string ma_name_full =
      "../out/mat/mat_" + maname + "_" + get_timestamp() + ".ma_dup";

#include "thinning.h"
  // compute face importances
  load_all_mat_face_importance_globally(mat, false /*is_sort_randomly*/,
                                        false /*is_debug*/);

  std::ofstream fout(ma_name_full);
  fout << mat.vertices->size() << " " << mat.edges.size() << " "
       << mat.faces.size() << " " << mat.tets.size() << std::endl;

  // we might have redundant vertices
  // but we need those indices, so store it as delete flag = true
  for (int i = 0; i < mat.vertices->size(); i++) {
    Vector3 pos = mat.vertices->at(i).center;
    fout << "v " << std::setiosflags(std::ios::fixed) << std::setprecision(15)
         << pos[0] << " " << pos[1] << " " << pos[2] << " "
         << mat.vertices->at(i).radius << " " << int(mat.vertices->at(i).type);

    // store this for not showing redundant vertices
    if (mat.vertices->at(i).is_deleted) {
      fout << " 1";
    } else {
      fout << " 0";
    }
    // store dup_cnt
    fout << " " << mat.vertices->at(i).dup_cnt;
    fout << std::endl;
  }

  for (int i = 0; i < mat.edges.size(); i++) {
    const auto& mat_e = mat.edges[i];
    // if (mat_e.is_deleted) continue;
    aint2 edge = {{mat_e.vertices_[0], mat_e.vertices_[1]}};
    std::sort(edge.begin(), edge.end());
    fout << "e " << edge[0] << " " << edge[1];
    // is this edge a feature edge
    if (mat_e.is_intf) {
      fout << " " << 1;  // internal
    } else if (mat_e.is_extf) {
      fout << " " << 2;  // external
    } else {
      fout << " " << 0;
    }
    fout << " " << mat_e.is_deleted;
    fout << " " << mat_e.dup_cnt;
    fout << std::endl;
  }

  for (int i = 0; i < mat.faces.size(); i++) {
    const auto& mat_f = mat.faces[i];
    // if (mat_f.is_deleted) continue;
    fout << "f";
    for (uint v = 0; v < 3; v++) fout << " " << mat_f.vertices_[v];
    fout << " " << std::setiosflags(std::ios::fixed) << std::setprecision(15)
         << mat_f.importance;
    fout << " " << mat_f.is_deleted;
    fout << " " << mat_f.dup_cnt;
    fout << std::endl;
  }

  for (int i = 0; i < mat.tets.size(); i++) {
    const auto& mat_t = mat.tets[i];
    fout << "t";
    for (int tet_vid : mat_t.vertices_) fout << " " << tet_vid;
    fout << " " << mat_t.is_deleted;
    fout << " " << mat_t.dup_cnt;
    fout << std::endl;
  }
  fout.close();

  printf("saved mat at: %s \n", ma_name_full.c_str());
}

/*
 * numVertices numEdges numFaces numTets
 * v x y z r flag_type flag_delete dup_cnt
 * e v1 v2 flag_type flag_delete dup_cnt
 * f v1 v2 v3 importance flag_delete dup_cnt
 * t v1 v2 v3 v4 flag_delete dup_cnt
 */
void load_ma_dup_cnt(const std::string& ma_path,
                     std::vector<MedialSphere>& all_medial_spheres,
                     MedialMesh& mat) {
  std::ifstream ma_file(ma_path);
  int num_vs, num_es, num_fs, num_ts;
  int type, is_deleted, dup_cnt;
  double importance;
  char ch;
  ma_file >> num_vs >> num_es >> num_fs >> num_ts;
  printf("loading .ma with num_vs: %d, num_es: %d, num_fs: %d, num_ts: %d\n",
         num_vs, num_es, num_fs, num_ts);

  // load spheres
  mat.clear();
  all_medial_spheres.clear();
  for (int i = 0; i < num_vs; i++) {
    MedialSphere msphere(all_medial_spheres.size(), Vector3(0, 0, 0),
                         Vector3(0, 0, 0), SphereType::T_2);
    ma_file >> ch >> msphere.center[0] >> msphere.center[1] >>
        msphere.center[2];
    ma_file >> msphere.radius;
    ma_file >> type;
    if (type != SphereType::T_UNK)
      msphere.type = SphereType(type);  // else as T2 sphere
    ma_file >> is_deleted;
    ma_file >> dup_cnt;
    msphere.is_deleted = is_deleted;
    msphere.dup_cnt = dup_cnt;
    all_medial_spheres.push_back(msphere);
    // mat
    if (!msphere.is_deleted) mat.numSpheres_active += dup_cnt;
  }
  mat.vertices = &all_medial_spheres;

  // load edge
  int e1, e2;
  for (int e = 0; e < num_es; e++) {
    ma_file >> ch >> e1 >> e2 >> type >> is_deleted >> dup_cnt;
    if (is_deleted) continue;
    // printf("ma e %d has (%d,%d)\n", e, e1, e2);
    int eid = mat.create_edge(e1, e2, dup_cnt);
    auto& medge = mat.edges.at(eid);
    if (type == 1) {
      medge.is_intf = true;
    } else if (type == 2) {
      medge.is_extf = true;
    }
    if (dup_cnt > 1) {
      printf("[MAT_DUP] edge %d (%d,%d) has dup_cnt %d\n", e, e1, e2, dup_cnt);
    }
  }
  // load facees
  int f1, f2, f3;
  for (int f = 0; f < num_fs; f++) {
    ma_file >> ch >> f1 >> f2 >> f3 >> importance >> is_deleted >> dup_cnt;
    // printf("ma f %d has (%d,%d,%d)\n", f, f1, f2, f3);
    aint3 fvs = {{f1, f2, f3}};
    int mfid = mat.create_face(fvs, dup_cnt);
    mat.faces.at(mfid).importance = importance;
    if (dup_cnt > 1) {
      printf("[MAT_DUP] face %d (%d,%d,%d) has dup_cnt %d\n", f, f1, f2, f3,
             dup_cnt);
    }
  }
  // load tets
  int t1, t2, t3, t4;
  for (int t = 0; t < num_ts; t++) {
    ma_file >> ch >> t1 >> t2 >> t3 >> t4 >> is_deleted >> dup_cnt;
    aint4 tvs = {{t1, t2, t3, t4}};
    mat.create_tet(tvs, dup_cnt);

    // error?
    if (dup_cnt > 1) {
      printf("[MAT_DUP] ERROR: tet %d (%d,%d,%d,%d) has dup_cnt %d\n", t, t1,
             t2, t3, t4, dup_cnt);
    }
  }
  ma_file.close();

  int euler = compute_Euler(mat);
  printf("[Euler] MedialMesh has Euler %d: v %ld, e %ld, f %ld, t %ld \n",
         euler, mat.vertices->size(), mat.edges.size(), mat.faces.size(),
         mat.tets.size());
}

// shibo: houdini output test (file .geo)
bool save_input_tet_houdini(const TetMesh& tet_mesh, std::string tet_name) {
  IO::Geometry geometry;
  IO::GeometryWriter geometry_writer(".");
  std::string tet_path = "../out/tet/tet_" + tet_name + "_" + get_timestamp();

  std::vector<float3> tet_vertices;
  for (int tvid = 0; tvid < tet_mesh.tet_vertices.size() / 3; tvid++) {
    float3 vertex = {tet_mesh.tet_vertices.at(tvid * 3),
                     tet_mesh.tet_vertices.at(tvid * 3 + 1),
                     tet_mesh.tet_vertices.at(tvid * 3 + 2)};
    tet_vertices.push_back(vertex);
  }

  std::vector<std::vector<unsigned>> tet_faces;
  std::vector<int> tet_faces_tid;
  int num_tets = tet_mesh.tet_indices.size() / 4;
  for (int tid = 0; tid < num_tets; tid++) {
    for (int j = 0; j < 4; j++) {
      int v1 = tet_mesh.tet_indices.at(tid * 4 + j);
      int v2 = tet_mesh.tet_indices.at(tid * 4 + (j + 1) % 4);
      int v3 = tet_mesh.tet_indices.at(tid * 4 + (j + 2) % 4);
      std::vector<unsigned> tmp = {v1, v2, v3};
      tet_faces.push_back(tmp);
      tet_faces_tid.push_back(tid);
    }
  }
  printf("tet_faces.size() %zu, num_tets: %d\n", tet_faces.size(), num_tets);

  assert(tet_faces.size() == num_tets * 4);
  assert(tet_faces_tid.size() == tet_faces.size());

  geometry.AddParticleAttribute("P", tet_vertices);
  geometry.AddPolygon(tet_faces);
  geometry.AddPrimitiveAttribute("PrimAttr", tet_faces_tid);
  printf("saving... \n");
  geometry_writer.OutputGeometry(tet_path, geometry);

  printf("saved tet houdini file: %s \n", tet_path.c_str());
}