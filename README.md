# MATTopo: Topology-preserving Medial Axis Transform with Restricted Power Diagram

- Project Link: https://ningnawang.github.io/projects/2024_mattopo
- Paper Link: https://arxiv.org/abs/2403.18761

## Please kindly cite our paper as:
```
@article{wang2024mattopo,
  title.    = {MATTopo: Topology-preserving Medial Axis Transform with Restricted Power Diagram},
  author.   = {Wang, Ningna and Huang, Hui and Song, Shibo and Wang, Bin and Wang, Wenping and Guo, Xiaohu},
  journal   = {ACM Transactions on Graphics (TOG)},
  year      = {2024},
  publisher = {ACM},
  booktitle = {ACM SIGGRAPH Asia 2024 Papers},
  series = {SIGGRAPH Asia '24}
}
```

## Dependencies:

1. Please make sure following libs are up-to-date:

```
apt install zlib1g-dev
sudo apt install mesa-common-dev
sudo apt-get install libcgal-dev
```

2. This project heavily relies on sub-modules in **extern/mat_modules**, please make sure **[mat_modules](https://github.com/ningnawang/mat_modules)** exists and is runnable. Normally you don't need to do anything with it, linux command `$cmake ..` will do the work for you (see Section **Build**).

3. Please also install following dependencies:
    * CMake
    * CGAL
    * Eigen3

## Build

This project has been tested on two Linux machines. 
Please build our code using **CMake**.

- Machine 1: GPU GeForce RTX 2080Ti
* Ubuntu 20.04
* Nvidia driver 525
* CUDA 11.8
* GCC 9.4.0
* CMake 3.16.3
* Geogram v1.7.5

- Machine 2: GPU GeForce RTX 4090
* Ubuntu 20.04
* Nvidia driver 535
* CUDA 12.2
* GCC 9.4.0
* CMake 3.16.3
* Geogram v1.7.5


For Linux, use the following steps to build:

```bash
mkdir build
cd build
cmake ..
make -j4
```

## Usage

The input of our method is the tetrahedral mesh in **.msh** format, one can use the [ftetwild](https://github.com/wildmeshing/fTetWild) to generate.

1. For running **NON-CAD** model (organic models), please run:
```bash
$ ./bin/MATTOPO <tet_mesh.msh> <nb_spheres> non_cad
```

2. For running **CAD** model with sharp edges and corners predetected, please run:
```bash
$ ./bin/MATTOPO <tet_mesh.msh> <nb_spheres>
```

## Output & Format

For the generated 3D medial mesh, there are three medial elements:
1. vertices, represents medial spheres.
2. edges, represents medial cones:
   - edges not shared by any faces (1D curves)
   - edges shared by at least one face
3. faces, represents medial slabs.


The result will be saved under the **out/mat/** folder. It saves the generated medial mesh in both **.ma** format and **.ply** format. 
We provide a GUI in [polyscope](https://polyscope.run/) which shows the generated medial mesh at the end of the program.

To open a saved 3D medial mesh in format `*.ma` or `*.ply`, we recommend following option:

1. The **.ma** file is the commonly used MAT format:
```
numVertices numEdges numFaces
v x y z r
e v1 v2
f v1 v2 v3
```
One can load the *.ma file using **Blender** with an open-sourced [blender-mat-addon](https://github.com/songshibo/blender-mat-addon).

2. The header ofr **.ply** will show the statistics of all vertices, faces and edges. Here edges consists of two parts as described above.
```
element vertex xxx
...
element face xxx
...
element edge xxx
```

We recommend using **Blender** to view the 3D medial meshes' in PLY format, as it shows correct statistics of not only vertices and faces, but also all edges, including 1d curves. If you open the PLY/OBJ file with **MeshLab**, only vertices and faces shows correctly, however, 1D curves will **NOT** be shown/counted in MeshLab.


