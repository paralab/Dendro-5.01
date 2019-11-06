/*
*@brief Contains utility functions for QUADGRAV simulation.
*/


#ifndef SFCSORTBENCH_GRUTILS_H
#define SFCSORTBENCH_GRUTILS_H

#include "point.h"
#include "parameters.h"
#include "mesh.h"
#include "block.h"
#include "parUtils.h"
#include "json.hpp"
#include "dendroProfileParams.h"
#include "profile_params.h"


#define Rx (quadgrav::QUADGRAV_COMPD_MAX[0]-quadgrav::QUADGRAV_COMPD_MIN[0])
#define Ry (quadgrav::QUADGRAV_COMPD_MAX[1]-quadgrav::QUADGRAV_COMPD_MIN[1])
#define Rz (quadgrav::QUADGRAV_COMPD_MAX[2]-quadgrav::QUADGRAV_COMPD_MIN[2])

#define RgX (quadgrav::QUADGRAV_OCTREE_MAX[0]-quadgrav::QUADGRAV_OCTREE_MIN[0])
#define RgY (quadgrav::QUADGRAV_OCTREE_MAX[1]-quadgrav::QUADGRAV_OCTREE_MIN[1])
#define RgZ (quadgrav::QUADGRAV_OCTREE_MAX[2]-quadgrav::QUADGRAV_OCTREE_MIN[2])

#define GRIDX_TO_X(xg) (((Rx/RgX)*(xg-quadgrav::QUADGRAV_OCTREE_MIN[0]))+quadgrav::QUADGRAV_COMPD_MIN[0])
#define GRIDY_TO_Y(yg) (((Ry/RgY)*(yg-quadgrav::QUADGRAV_OCTREE_MIN[1]))+quadgrav::QUADGRAV_COMPD_MIN[1])
#define GRIDZ_TO_Z(zg) (((Rz/RgZ)*(zg-quadgrav::QUADGRAV_OCTREE_MIN[2]))+quadgrav::QUADGRAV_COMPD_MIN[2])

#define X_TO_GRIDX(xc) (((RgX/Rx)*(xc-quadgrav::QUADGRAV_COMPD_MIN[0]))+quadgrav::QUADGRAV_OCTREE_MIN[0])
#define Y_TO_GRIDY(yc) (((RgY/Ry)*(yc-quadgrav::QUADGRAV_COMPD_MIN[1]))+quadgrav::QUADGRAV_OCTREE_MIN[1])
#define Z_TO_GRIDZ(zc) (((RgZ/Rz)*(zc-quadgrav::QUADGRAV_COMPD_MIN[2]))+quadgrav::QUADGRAV_OCTREE_MIN[2])


using json = nlohmann::json;
namespace quadgrav
{
/**
 * @brief These variable indexes are based on the variables defined in rkQUADGRAV.h
 * */
enum VAR {U_CHI=0,U_PHI};

static const char * QUADGRAV_VAR_NAMES[]={"U_CHI","U_PHI"};

/**
 * @brief internal variables needed for rk update.
 * */


 /**
  * @brief: Read the parameter file and initialize the variables in parameters.h file.
  * @param[in] fName: file name
  * @param[in] comm: MPI communicator.
  * */
  void readParamFile(const char * fName,MPI_Comm comm);


/**
 * @brief Initialize all the variables for a given point in space.
 * @param [in] coord: coordinates of the point.
 * @param [out] var: pointer to the list of variables, computed. var size should be (VAR::U_SYMAT5+1)
 * @note This function is taken from the old single core quadgrav version.
 **/

 // Initial data
 void initData(const double xx1,const double yy1,const double zz1, double *var);

 /**
  * analytical solution based on the d'Alembert's formular.
  *
  * */
 void analyticalSol(const double xx1,const double yy1,const double zz1,const double t, double *var);

 /**
  * @brief: Generates block adaptive octree for the given binary blockhole problem.
  *
  * */

  void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,const Point& pt_min,const Point & pt_max,const unsigned int regLev,const unsigned int maxDepth,MPI_Comm comm);

  /**
   * @brief wavelet tolerance as a function of space.
   * */
  double computeWTol(double x,double y,double z,double tol_min);

}// end of namespace quadgrav



namespace quadgrav
{

    namespace timer
    {

        /**@brief initialize all the flop counters. */
        void initFlops();

        /**@brief clears the snapshot counter for time profiler variables*/
        void resetSnapshot();


       /**@brief reduce min mean max.
        * @param [in] stat: local time
        * @param [out] stat_g 0-min, 1-mean 2-max
       * */
       template<typename T>
       void computeOverallStats(T *stat, T *stat_g, MPI_Comm comm)
       {
           int rank,npes;
           MPI_Comm_size(comm,&npes);
           MPI_Comm_rank(comm,&rank);

           par::Mpi_Reduce(stat,stat_g,1,MPI_MIN,0,comm);
           par::Mpi_Reduce(stat,stat_g+1,1,MPI_SUM,0,comm);
           par::Mpi_Reduce(stat,stat_g+2,1,MPI_MAX,0,comm);
           stat_g[1]/=(npes);

       }


        /** @breif : printout the profile parameters. */
        void profileInfo(const char* filePrefix,const ot::Mesh* pMesh);

        /** @breif : printout the profile parameters (intermediate profile information). */
        void profileInfoIntermediate(const char* filePrefix,const ot::Mesh* pMesh,const unsigned int currentStep);


    }


}

#endif //SFCSORTBENCH_GRUTILS_H
