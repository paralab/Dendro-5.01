//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief This file contains all the parameters related to QUADGRAV simulation.
*/
//

#ifndef SFCSORTBENCH_PARAMETERS_H
#define SFCSORTBENCH_PARAMETERS_H

#include <string.h>
#include <iostream>




namespace quadgrav
{

    /**@brief element order*/
    static const unsigned int QUADGRAV_ELE_ORDER=4;

    /**@brief number of variables*/
    static const unsigned int QUADGRAV_NUM_VARS=2;

    /***@brief number of RK45 stages*/
    static const unsigned int QUADGRAV_RK45_STAGES=6;

    /***@brief number of RK4 stages*/
    static const unsigned int QUADGRAV_RK4_STAGES=4;

    /**@brief number of rk4 stages*/
    static const unsigned int QUADGRAV_RK3_STAGES=3;
    
    /**@brief: parameter used for adaptive time step update. */
    static const double QUADGRAV_SAFETY_FAC=0.8;

    /**@brief number of internal variables*/
    static const unsigned int QUADGRAV_NUM_VARS_INTENL=(QUADGRAV_RK45_STAGES+1)*QUADGRAV_NUM_VARS;

    /**@brief min bh domain add these to the parameter file.*/
    extern double QUADGRAV_COMPD_MIN[3];
    /**@brief min bh domain @todo add these to the parameter file. */
    extern double QUADGRAV_COMPD_MAX[3];

    /**@brief CFL stability number number (specifies how dt=QUADGRAV_CFL_FACTOR*dx)*/
    extern double QUADGRAV_CFL_FACTOR;


    /**@brief min coords of the OCTREE */
    extern double QUADGRAV_OCTREE_MIN[3];
    /**@brief max coords of the OCTREE */
    extern double QUADGRAV_OCTREE_MAX[3];

    /**@brief solution output frequency*/
    extern unsigned int QUADGRAV_IO_OUTPUT_FREQ;

    /**@brief timestep norms out put freq.*/
    extern unsigned int QUADGRAV_TIME_STEP_OUTPUT_FREQ;

    /**@brief remesh test frequency*/
    extern unsigned int QUADGRAV_REMESH_TEST_FREQ;

    /**@brief checkpoint store frequency*/
    extern unsigned int QUADGRAV_CHECKPT_FREQ;

    /**@brief restore the solver from check point if set to 1. */
    extern unsigned int QUADGRAV_RESTORE_SOLVER;

    /**@brief use the block adaptivity and disable the AMR*/
    extern unsigned int QUADGRAV_ENABLE_BLOCK_ADAPTIVITY;

    /**@brief file prefix for VTU*/
    extern std::string QUADGRAV_VTU_FILE_PREFIX;

    /**@brief file prefix for write check point*/
    extern std::string QUADGRAV_CHKPT_FILE_PREFIX;

    /**@brief file prefix to write profile info.*/
    extern std::string QUADGRAV_PROFILE_FILE_PREFIX;

    /**@brief number of refine variables*/
    extern unsigned int QUADGRAV_NUM_REFINE_VARS;

    /**@brief indices of refine var ids*/
    extern unsigned int QUADGRAV_REFINE_VARIABLE_INDICES[QUADGRAV_NUM_VARS];

    /**@brief number of evolution variables written to vtu files*/
    extern unsigned int QUADGRAV_NUM_EVOL_VARS_VTU_OUTPUT;

    /**@brief evolution variable IDs written to vtu files*/
    extern unsigned int QUADGRAV_VTU_OUTPUT_EVOL_INDICES[QUADGRAV_NUM_VARS];

    /**@brief solution output gap (instead of freq. we can use to output the solution if currentTime > lastIOOutputTime + QUADGRAV_IO_OUTPUT_GAP)*/
    extern  double QUADGRAV_IO_OUTPUT_GAP;

    /**@brief prefered grain sz to use when selecting active npes*/
    extern unsigned int QUADGRAV_DENDRO_GRAIN_SZ;

    /**@brief AMR coarsening factor (we coarsen if tol<QUADGRAV_DENDRO_AMR_FAC*QUADGRAV_WAVELET_TOL)*/
    extern double QUADGRAV_DENDRO_AMR_FAC;

    /**@brief wavelet tolerance value. */
    extern  double QUADGRAV_WAVELET_TOL;
    /**@brief load-imbalance tolerance value. */
    extern  double QUADGRAV_LOAD_IMB_TOL;
    /**@brief: Splitter fix value*/
    extern unsigned int QUADGRAV_SPLIT_FIX;

    /**@brief: async. communication at a time. (upper bound shoud be QUADGRAV_NUM_VARS) */
    extern unsigned int QUADGRAV_ASYNC_COMM_K;


    /**@brief simulation begin time. */
    extern double QUADGRAV_RK45_TIME_BEGIN;
    /**@brief simulation end time*/
    extern double QUADGRAV_RK45_TIME_END;
    /**@brief rk time step size. */
    extern double QUADGRAV_RK45_TIME_STEP_SIZE;

    /** desired tolerance value for the rk45 method (adaptive time stepping. )*/
    extern double QUADGRAV_RK45_DESIRED_TOL;

    /**@brief BBH initial data type */
    extern unsigned int QUADGRAV_ID_TYPE;

    /**@brief physical coordinates for grid, x_min */
    extern double QUADGRAV_GRID_MIN_X;

    /**@brief physical coordinates for grid, x_max */
    extern double QUADGRAV_GRID_MAX_X;

    /**@brief physical coordinates for grid, y_min */
    extern double QUADGRAV_GRID_MIN_Y;

    /**@brief physical coordinates for grid, y_max */
    extern double QUADGRAV_GRID_MAX_Y;

    /**@brief physical coordinates for grid, z_min */
    extern double QUADGRAV_GRID_MIN_Z;

    /**@brief physical coordinates for grid, z_max */
    extern double QUADGRAV_GRID_MAX_Z;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double QUADGRAV_BLK_MIN_X;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double QUADGRAV_BLK_MIN_Y;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double QUADGRAV_BLK_MIN_Z;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double QUADGRAV_BLK_MAX_X;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double QUADGRAV_BLK_MAX_Y;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double QUADGRAV_BLK_MAX_Z;

    /**@brief: dimension of the grid*/
    extern unsigned int QUADGRAV_DIM;

    /**@brief: max refinement level*/
    extern unsigned int QUADGRAV_MAXDEPTH;

    /**@brief: Kreiss-Oliger dissipation */
    extern double KO_DISS_SIGMA;

    /**@brief: Kreiss-Oliger dissipation */
    extern double KO_DISS_SIGMA;



    /**@brief: Initial data Gaussian amplitude */
    extern double QUADGRAV_ID_AMP1;

    /**@brief: Initial data Gaussian amplitude */
    extern double QUADGRAV_ID_AMP2;

    /**@brief: Initial data Gaussian width */
    extern double QUADGRAV_ID_DELTA1;

    /**@brief: Initial data Gaussian width */
    extern double QUADGRAV_ID_DELTA2;

    /**@brief: Initial data Gaussian x offset */
    extern double QUADGRAV_ID_XC1;
    extern double QUADGRAV_ID_YC1;
    extern double QUADGRAV_ID_ZC1;

    /**@brief: Initial data Gaussian x offset */
    extern double QUADGRAV_ID_XC2;
    extern double QUADGRAV_ID_YC2;
    extern double QUADGRAV_ID_ZC2;

    /**@brief: Initial data Gaussian elliptic x factor */
    extern double QUADGRAV_ID_EPSX1;

    /**@brief: Initial data Gaussian elliptic y factor */
    extern double QUADGRAV_ID_EPSY1;

    /**@brief: Initial data Gaussian elliptic z factor */
    extern double QUADGRAV_ID_EPSZ1;

    /**@brief: Initial data Gaussian elliptic x factor */
    extern double QUADGRAV_ID_EPSX2;

    /**@brief: Initial data Gaussian elliptic y factor */
    extern double QUADGRAV_ID_EPSY2;

    /**@brief: Initial data Gaussian elliptic z factor */
    extern double QUADGRAV_ID_EPSZ2;

    /**@brief: Initial data Gaussian R */
    extern double QUADGRAV_ID_R1;

    /**@brief: Initial data Gaussian R */
    extern double QUADGRAV_ID_R2;

    /**@brief: Initial data Gaussian nu */
    extern double QUADGRAV_ID_NU1;

    /**@brief: Initial data Gaussian nu */
    extern double QUADGRAV_ID_NU2;

    /**@brief: Initial data Gaussian Omega */
    extern double QUADGRAV_ID_OMEGA;

    /**@brief: wave speed direction x*/
    extern double QUADGRAV_WAVE_SPEED_X;

    /**@brief: wave speed direction y*/
    extern double QUADGRAV_WAVE_SPEED_Y;
    
    /**@brief: wave speed direction z*/
    extern double QUADGRAV_WAVE_SPEED_Z;
}

#endif //SFCSORTBENCH_PARAMETERS_H
