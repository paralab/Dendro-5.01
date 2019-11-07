# Quadratic Gravity project in Dendro

Hyun Lim

In this project, we study black hole dynamics in 
quadractic gravity (QG). All detailed discussion about
QG can be found [here](https://github.com/hlim88/quadGrav)

## Detailed instruction Build Dendro 

Dendro can be installed anywhere in your system but we suggest
that all repositories are downloaded in DENDRO root directory 
such as `${HOME}/DENDRO`

### Suggested directory structure

We recommend to have following directory structure:

```{engine=sh}
  ${HOME}/DENDRO
  ├── Dendro-5.01
  │   └── build
  └── local
      ├── bin
      ├── include
      ├── lib
      ├── lib64
      └── share
```

Below we use `${HOME}/DENDRO/local` for an installation directory.
Make sure to set your CMAKE prefix to this location:

    % export CMAKE_PREFIX_PATH=${HOME}/DENDRO/local

### Prerequisites

You will need the following tools:

- C++11 - capable compiler, such as gcc version >= 4.8;
- MPI libraries, compiled with the gcc compiler above and multithread support
  (`--enable-mpi-thread-multiple` for OpenMPI and
   `--enable-threads=multiple` for MPICH);
- cmake version > 2.8;
- Python version > 3.5; (For RHS generation scripts. Not required for 
building code unless you want to generate RHS during runtime)
- zlib compression library;
- Optional : BLAS and LAPACK

#### Building Dendro

Clone the master branch from the Dendro-5.01 git repo:
```{engine=sh}
   cd $HOME/DENDRO
   git clone https://github.com/paralab/Dendro-5.01.git
```    

After cloning the repo, follow the configure command:
(copy and paste below into your temrinal)
```{engine=sh}
   # in ${HOME}/DENDRO/Dendro-5.01:
   mkdir build ; cd build
   export CMAKE_PREFIX_PATH=${HOME}/DENDRO/local
   cmake .. \
       -DCMAKE_INSTALL_PREFIX=$CMAKE_PREFIX_PATH \
       -DALLTOALLV_FIX=OFF                       \
       -DHILBERT_ORDERING=ON                     \
       -DPROFILE_TREE_SORT=OFF                   \
       -DSPLITTER_SELECTION_FIX=OFF              \
       -DNUM_NPES_THRESHOLD=2                    \
```

Build and install inside of your `build` directory by:
```{engine=sh}
     make -j 
     make install
```
Note that if your machine do not support parallel job, 
you may not use `-j` option during `make`

### Building Dendro with Spack

This will be added once Dendro moves into Spack package repo.

### Running Dendro application : Concentrating with QG project

Once you successfully building the Dendro, you will have `quadgravSolver` 
in your `build` directory and/or `${HOME}/DENDRO/local/bin`.

You can run `quadgravSolver` as follow
```{english=sh}
  mpirun -np <number of mpi tasks> ./quadgravSolver <parameter file name>.par
```
An example parameter file can be found in `QuadGrav/par`. You can 
create your own parameter file based on this
