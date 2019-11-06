//
// Created by milinda on 8/23/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "parameters.h"

namespace quadgrav
{
    unsigned int QUADGRAV_IO_OUTPUT_FREQ=10;
    unsigned int QUADGRAV_TIME_STEP_OUTPUT_FREQ=10;
    unsigned int QUADGRAV_REMESH_TEST_FREQ=10;
    double QUADGRAV_IO_OUTPUT_GAP=1.0;

    double QUADGRAV_WAVELET_TOL=0.0001;

    double QUADGRAV_LOAD_IMB_TOL=0.1;
    unsigned int QUADGRAV_SPLIT_FIX=2;
    unsigned int QUADGRAV_ASYNC_COMM_K=4;
    double QUADGRAV_RK45_TIME_BEGIN=0;
    double QUADGRAV_RK45_TIME_END=10;

    double QUADGRAV_RK45_DESIRED_TOL=1e-6;


    unsigned int QUADGRAV_CHECKPT_FREQ=10;
    unsigned int QUADGRAV_RESTORE_SOLVER=0;
    unsigned int QUADGRAV_ENABLE_BLOCK_ADAPTIVITY=0;

    std::string QUADGRAV_VTU_FILE_PREFIX="quadgrav_gr";
    std::string QUADGRAV_CHKPT_FILE_PREFIX="quadgrav_cp";
    std::string QUADGRAV_PROFILE_FILE_PREFIX="quadgrav_prof";


    unsigned int QUADGRAV_DIM=3;
    unsigned int QUADGRAV_MAXDEPTH=8;

    unsigned int QUADGRAV_ID_TYPE=0;
    double QUADGRAV_CFL_FACTOR=0.1;

    double QUADGRAV_GRID_MIN_X=-50.0;
    double QUADGRAV_GRID_MAX_X=50.0;
    double QUADGRAV_GRID_MIN_Y=-50.0;
    double QUADGRAV_GRID_MAX_Y=50.0;
    double QUADGRAV_GRID_MIN_Z=-50.0;
    double QUADGRAV_GRID_MAX_Z=50.0;

    double QUADGRAV_BLK_MIN_X=-6.0;
    double QUADGRAV_BLK_MIN_Y=-6.0;
    double QUADGRAV_BLK_MIN_Z=-6.0;

    double QUADGRAV_BLK_MAX_X=6.0;
    double QUADGRAV_BLK_MAX_Y=6.0;
    double QUADGRAV_BLK_MAX_Z=6.0;

    double QUADGRAV_COMPD_MIN[3]={QUADGRAV_GRID_MIN_X,QUADGRAV_GRID_MIN_Y,QUADGRAV_GRID_MIN_Z};
    double QUADGRAV_COMPD_MAX[3]={QUADGRAV_GRID_MAX_X,QUADGRAV_GRID_MAX_Y,QUADGRAV_GRID_MAX_Z};

    double QUADGRAV_OCTREE_MIN[3]={0.0,0.0,0.0};
    double QUADGRAV_OCTREE_MAX[3]={(double)(1u<<QUADGRAV_MAXDEPTH),(double)(1u<<QUADGRAV_MAXDEPTH),(double)(1u<<QUADGRAV_MAXDEPTH)};

    //@note assumes the computational domain is a cube as well.
    double QUADGRAV_RK45_TIME_STEP_SIZE=QUADGRAV_CFL_FACTOR*(QUADGRAV_COMPD_MAX[0]-QUADGRAV_COMPD_MIN[0])*(1.0/(double)(1u<<QUADGRAV_MAXDEPTH));

    double KO_DISS_SIGMA=0.01;


    unsigned int QUADGRAV_DENDRO_GRAIN_SZ=1000;

    double QUADGRAV_DENDRO_AMR_FAC=0.1;

    unsigned int QUADGRAV_NUM_REFINE_VARS=2;
    unsigned int QUADGRAV_REFINE_VARIABLE_INDICES[QUADGRAV_NUM_VARS]={0,1};

    unsigned int QUADGRAV_NUM_EVOL_VARS_VTU_OUTPUT=2;
    unsigned int QUADGRAV_VTU_OUTPUT_EVOL_INDICES[QUADGRAV_NUM_VARS]={0,1};

    double QUADGRAV_ID_AMP1 = 0.5;
    double QUADGRAV_ID_AMP2 = 0.5;
    double QUADGRAV_ID_DELTA1 = 1.0;
    double QUADGRAV_ID_DELTA2 = 1.0;
    double QUADGRAV_ID_XC1 = 0.0;
    double QUADGRAV_ID_YC1 = 0.0;
    double QUADGRAV_ID_ZC1 = 0.0;
    double QUADGRAV_ID_XC2 = 0.0;
    double QUADGRAV_ID_YC2 = 0.0;
    double QUADGRAV_ID_ZC2 = 0.0;
    double QUADGRAV_ID_EPSX1 = 1.0;
    double QUADGRAV_ID_EPSY1 = 1.0;
    double QUADGRAV_ID_EPSZ1 = 1.0;
    double QUADGRAV_ID_EPSX2 = 1.0;
    double QUADGRAV_ID_EPSY2 = 1.0;
    double QUADGRAV_ID_EPSZ2 = 1.0;
    double QUADGRAV_ID_R1 = 0.0;
    double QUADGRAV_ID_R2 = 0.0;
    double QUADGRAV_ID_NU1 = 0.0;
    double QUADGRAV_ID_NU2 = 0.0;
    double QUADGRAV_ID_OMEGA = 0.0;

    double QUADGRAV_WAVE_SPEED_X = 1.0;
    double QUADGRAV_WAVE_SPEED_Y = 0.0;
    double QUADGRAV_WAVE_SPEED_Z = 0.0;
    
}
