<img src="https://github.com/paralab/Dendro-5.01/blob/master/docs/fig/dendro.png" alt="Dendro" width="400"/>

## What is Dendro ?

"Dendro" in Greek language means tree. The Dendro library is a large scale (262K cores on ORNL's Titan) distributed memory 
adaptive octree framework. The main goal of Dendro is to perform large scale multiphysics simulations efficeiently in mordern supercomputers. Dendro consists of efficient parallel data structures and algorithms to perform variational ( finite element) methods and finite difference mthods on 2:1 balanced arbitary adaptive octrees which enables the users to perform simulations raning from black holes (binary black hole mergers) to blood flow in human body, where applications ranging from relativity, astrophysics to biomedical engineering.  

***

## Get Dendro

You can clone the repository using , `git clone https://github.com/paralab/Dendro-5.01.git`

## How to build Dendro-5.0 ?

You need CMake to build dendro. Create a build directory using 'mkdir build'. Then go into the build directory by 'cd build' then execute 'ccmake ..' to generate the make files. You can build Dendro-5.0 with several options. 

* `ALLTOALLV_FIX` : OFF,Need to turn off
* `DIM_2`: OFF, This can be turned on if you need to run Dendro-5.0 in 2D case. default: OFF (Which means it assumes 3D domain) 
* `HILBERT_ORDERING`:ON, This specify which SFC to use to partition the data. HILBERT_ORDERING: ON means it uses Hilbert curve, otherwise it uses Morton curve for partitioning. 
* `PROFILE_TREE_SORT`: OFF 
* `NUM_NPES_THRESHOLD`: square root of P (number of processors) 
* `SPLITTER_SELECTION_FIX`: ON. This will perform the data exchange in the octree partitioning in stages. This is mandatory when you run dendro in very large scale. 

## Simple simulation: Nonlinear Sigma Model (NLSigma)

NlSigma folder consists of simple, non lineat wave equation with adaptive mesh refinement (AMR). You can copy the parameter file from `NLSigma/pars` folder and simply run `mpirun -np 8 ./NlSigma/nlsmSolver nlsm.par.json`, on  your lattop to large supercomputer with higher resolution. 

<div id="mainDiv">
 <div id="divTwo" class="boxes">
	 <img src="https://github.com/paralab/Dendro-5.01/blob/master/docs/fig/nlsmB7.png" alt="nlsm" width="200"/>
    </div>
 <div id="divTwo" class="boxes">
	<img src="https://github.com/paralab/Dendro-5.01/blob/master/docs/fig/nlsmB11.png" alt="nlsm" width="200"/>
    </div>
 <div id="divTwo" class="boxes">
	<img src="https://github.com/paralab/Dendro-5.01/blob/master/docs/fig/nlsmB16.png" alt="nlsm" width="200"/>
    </div>
 <div id="divTwo" class="boxes">
	<img src="https://github.com/paralab/Dendro-5.01/blob/master/docs/fig/nlsmB44.png" alt="nlsm" width="200"/>
    </div>
</div>

You can write the equations in symbolic python which generate the C compute kernel. Look at `nlsm.py` 

```
import dendro
from sympy import *
###############################################################
#  initialize
###############################################################
r = symbols('r')
# declare functions
chi = dendro.scalar("chi","[pp]")
phi = dendro.scalar("phi","[pp]")
d = dendro.set_first_derivative('grad')    # first argument is direction
d2s = dendro.set_second_derivative('grad2')  # first 2 arguments are directions
d2 = dendro.d2

###############################################################
#  evolution equations
###############################################################

phi_rhs = sum( d2(i,i,chi)  for i in dendro.e_i ) - sin(2*chi)/r**2
chi_rhs = phi

###############################################################
#  evolution equations
###############################################################
outs = [phi_rhs, chi_rhs]
vnames = ['phi_rhs', 'chi_rhs']
dendro.generate(outs, vnames, '[pp]')
```

Which generate the code, 

```
// Dendro: {{{ 
// Dendro: original ops:  10
// Dendro: printing temp variables

// Dendro: printing variables
//--
phi_rhs[pp] = grad2_0_0_chi[pp] + grad2_1_1_chi[pp] + grad2_2_2_chi[pp] - sin(2*chi[pp])/pow(r, 2);
//--
chi_rhs[pp] = phi[pp];
// Dendro: reduced ops:  10
// Dendro: }}}
```


***

## Scalability on octree generation and partitioning. 

We have performed octree generation and partitioning up to 262144 cores in ORNL's titan super computer. We have managed to partition 1.3x10^12 octants among 262144 processors with in 4 seconds.

<img src="https://github.com/paralab/Dendro-5.01/blob/master/docs/fig/titan_ws.png" alt="weak scaling on octrees" width="800"/>
<img src="https://github.com/paralab/Dendro-5.01/blob/master/docs/fig/titan_ss.png" alt="strong scaling on octrees" width="800"/>


## Publications
* Milinda Fernando, David Neilsen, Hyun Lim, Eric Hirschmann, Hari Sundar, ”Massively Parallel Simulations of Binary Black Hole Intermediate-Mass-Ratio Inspirals” SIAM Journal on Scientific Computing 2019. 'https://doi.org/10.1137/18M1196972'
* Milinda Fernando, David Neilsen, Hari Sundar, ”A scalable framework for Adaptive Computational General Relativity on Heterogeneous Clusters”, (ACM International Conference on Supercomputing, ICS’19)
* Milinda Fernando, Dmitry Duplyakin, and Hari Sundar. 2017. ”Machine and Application Aware Partitioning for Adaptive Mesh Refinement Applications”. In Proceedings of the 26th International Symposium on High-Performance Parallel and Distributed Computing (HPDC ’17). ACM, New York, NY, USA, 231-242. DOI: 'https://doi.org/10.1145/3078597.3078610'
 
