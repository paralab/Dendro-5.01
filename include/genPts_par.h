/*
 *
 * @author: Milinda Fernando
 * School of Computing, University of Utah
 * @date: 11/10/2015
 *
 *
 * Contains the Normal(Gaussian) and Logarithmic normal random number (octants)
 * generator based on new c+11 rand engines.
 *
 *
 * */

#ifndef GEN_GAUSS_H
#define GEN_GAUSS_H

#include <mpi.h>

#include <chrono>
#include <iostream>
#include <random>

#include "TreeNode.h"
#include "dendro.h"

void genGauss(const double& sd, const DendroIntL numPts, int dim,
              char* filePrefix, MPI_Comm comm);
void genGauss(const double& sd, const DendroIntL numPts, int dim, double* xyz);
void genLogarithmicGauss(const double& sd, const DendroIntL numPts, int dim,
                         char* filePrefix);
void genUniformRealDis(const double& sd, const DendroIntL numPts, int dim,
                       double* xyz);
void genLogarithmicGauss(const double& sd, const DendroIntL numPts, int dim,
                         double* xyz);
void pts2Octants(std::vector<ot::TreeNode>& pNodes, double* pts,
                 DendroIntL totPts, unsigned int dim, unsigned int maxDepth);

#endif
