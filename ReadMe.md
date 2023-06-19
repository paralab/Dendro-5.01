For questions: [![Gitter](https://badges.gitter.im/Dendro-5-01/community.svg)](https://gitter.im/Dendro-5-01/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

[![DOI](https://zenodo.org/badge/102657919.svg)](https://zenodo.org/badge/latestdoi/102657919)

<img src="https://github.com/paralab/Dendro-5.01/blob/master/docs/fig/dendro.png" alt="Dendro" width="400"/>

## What is Dendro ?

"Dendro" in Greek language means tree. The Dendro library is a large scale (262K cores on ORNL's Titan) distributed memory
adaptive octree framework. The main goal of Dendro is to perform large scale multiphysics simulations efficiently in mordern supercomputers. Dendro consists of efficient parallel data structures and algorithms to perform variational (finite element) methods and finite difference mthods on 2:1 balanced arbitary adaptive octrees which enables the users to perform simulations raning from black holes (binary black hole mergers) to blood flow in human body, where applications ranging from relativity, astrophysics to biomedical engineering.

---

## This Repository

This repository contains the C++ Dendro library. There are a few examples included, but the majority of the simulations that use this library are housed in a separate repository. You can check out [Dendro-GR](https://github.com/paralab/Dendro-GR), a large portion of the General Relativity project that does the full physics simulations. In particular, be sure to check out the examples in the README, such as the `NLSigma` simulation.

## Get Dendro

You can clone the repository using `git clone https://github.com/paralab/Dendro-5.01.git`

## How do you build Dendro-5.0?

You need CMake to build dendro. Create a build directory using `mkdir build`. Then go into the build directory by `cd build` then execute `ccmake ..` to generate the make files. You can build Dendro-5.0 with several options (via `ccmake` or CMake flags of the same name).

- `ALLTOALLV_FIX`: `OFF`, This must be set to `OFF`.
- `DIM_2`: `OFF`, This can be turned on if you need to run Dendro-5.0 in a 2D simulation. When set to `OFF`, it assumes a 3D domain.
- `HILBERT_ORDERING`: `ON`, This specify which SFC to use to partition the data. Seetting it to `ON` means it will use the Hilbert curve, otherwise it will use the Morton curve for partitioning.
- `PROFILE_TREE_SORT`: `OFF`
- `NUM_NPES_THRESHOLD`: the square root of `P` (number of processors)
- `SPLITTER_SELECTION_FIX`: `ON`. This will perform the data exchange in the octree partitioning in stages. This is mandatory when you run Dendro in very large scale.

## Using Dendro in a project

We are in the process of writing more complete documentation for how to build a project with Dendro. This repository should be included as a whole submodule or as its own folder.

You can see the [Dendro-GR](https://github.com/paralab/Dendro-GR) repository (especially the `NLSigma` case) for multiple working examples of how Dendro is incorporated.

---

## Scalability on octree generation and partitioning.

We have performed octree generation and partitioning up to 262144 cores in ORNL's titan super computer. We have managed to partition 1.3x10^12 octants among 262144 processors with in 4 seconds.

<img src="https://github.com/paralab/Dendro-5.01/blob/master/docs/fig/titan_ws.png" alt="weak scaling on octrees" width="800"/>
<img src="https://github.com/paralab/Dendro-5.01/blob/master/docs/fig/titan_ss.png" alt="strong scaling on octrees" width="800"/>

## Publications

- Milinda Fernando, David Neilsen, Hyun Lim, Eric Hirschmann, Hari Sundar, ”Massively Parallel Simulations of Binary Black Hole Intermediate-Mass-Ratio Inspirals” SIAM Journal on Scientific Computing 2019. 'https://doi.org/10.1137/18M1196972'
- Milinda Fernando, David Neilsen, Hari Sundar, ”A scalable framework for Adaptive Computational General Relativity on Heterogeneous Clusters”, (ACM International Conference on Supercomputing, ICS’19)
- Milinda Fernando, Dmitry Duplyakin, and Hari Sundar. 2017. ”Machine and Application Aware Partitioning for Adaptive Mesh Refinement Applications”. In Proceedings of the 26th International Symposium on High-Performance Parallel and Distributed Computing (HPDC ’17). ACM, New York, NY, USA, 231-242. DOI: 'https://doi.org/10.1145/3078597.3078610'
