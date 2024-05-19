#include <polyscope/polyscope.h>
#include <polyscope/volume_mesh.h>
#include <stdio.h>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "../extern/mshloader/MshLoader.h"

using namespace PyMesh;
bool read_tet_mesh(const std::string& filename,
                   std::vector<std::array<float, 3>>& voro_points,
                   std::vector<std::array<int, 4>>& voro_tets,
                   std::vector<int>& voro_tets_sites) {
  int s;
  int n_pts, n_tets;
  std::string ext = filename.substr(filename.find_last_of('.') + 1);
  // printf("loading file %s with ext %s", filename.c_str(), ext.c_str());
  std::cout << "loading file " << filename.c_str() << " with ext "
            << ext.c_str() << std::endl;
  if (ext == "msh") {
    MshLoader msh_loader(filename);
    std::vector<float> vs = msh_loader.get_nodes();
    std::vector<int> tets = msh_loader.get_elements();
    n_pts = vs.size() / 3;
    n_tets = tets.size() / 4;
    voro_points.resize(n_pts);
    voro_tets.resize(n_tets);
    voro_tets_sites.resize(n_tets);

    for (uint i = 0; i < n_pts; ++i) {
      for (uint j = 0; j < 3; j++) voro_points[i][j] = vs[i * 3 + j];
      // std::cout << "vertex " << i << ": (" << voro_points[i][0] << ", "
      //           << voro_points[i][1] << ", " << voro_points[i][2] << ")"
      //           << std::endl;
      // std::cout << "vertex " << i << ": (" << vs[i * 3 + 0] << ", "
      //           << vs[i * 3 + 1] << ", " << vs[i * 3 + 2] << ")" <<
      //           std::endl;
    }

    for (uint i = 0; i < n_tets; ++i) {
      for (uint j = 0; j < 4; j++) voro_tets[i][j] = tets[i * 4 + j];
      // std::cout << "tet " << i << ": (" << voro_tets[i][0] << ", "
      //           << voro_tets[i][1] << ", " << voro_tets[i][2] << ", "
      //           << voro_tets[i][3] << ")" << std::endl;
    }

  } else if (ext == "tet") {
    std::ifstream input(filename);
    if (input.fail()) return false;
    input >> n_pts >> n_tets;
    voro_points.resize(n_pts);
    voro_tets.resize(n_tets);
    voro_tets_sites.resize(n_tets, 0);

    for (uint i = 0; i < n_pts; ++i)
      input >> voro_points[i][0] >> voro_points[i][1] >> voro_points[i][2];

    for (uint i = 0; i < n_tets; ++i) {
      input >> s >> voro_tets[i][0] >> voro_tets[i][1] >> voro_tets[i][2] >>
          voro_tets[i][3];

      // when a tet is given a voronoi site
      if (s == 5) input >> voro_tets_sites[i];
    }
    input.close();
  } else {
    return false;
  }

  std::cout << "loaded voro_points: " << voro_points.size()
            << ", voro_tets: " << voro_tets.size() << std::endl;

  return true;
}

int main(int argc, char** argv) {
  // a few camera options
  polyscope::view::upDir = polyscope::UpDir::ZUp;

  // Initialize Polyscope
  polyscope::init();

  // Read mesh from file
  std::vector<std::array<float, 3>> voro_points;
  std::vector<std::array<int, 4>> voro_tets;
  std::vector<int> voro_tets_sites;

  if (!read_tet_mesh(argv[1], voro_points, voro_tets, voro_tets_sites)) {
    std::cerr << argv[0] << ": could not load file at " << argv[1] << std::endl;
    return 1;
  }

  // Register the volume mesh with Polyscope
  polyscope::registerTetMesh("my mesh", voro_points, voro_tets);

  // Add a scalar quantity
  polyscope::getVolumeMesh("my mesh")->addCellScalarQuantity("site",
                                                             voro_tets_sites);

  // Show the GUI
  polyscope::show();
}