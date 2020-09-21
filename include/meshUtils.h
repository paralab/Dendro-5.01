/**
 * @file meshUtils.h
 * @author Milinda Fernando
 * @brief Contains utility function related to easily generate a mesh. 
 * @version 0.1
 * @date 2020-01-01
 * 
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2020
 * 
 */


#pragma once
#include "mesh.h"
#include "octUtils.h"
#include <vector>
#include <iostream>
#include <functional>
#include "parUtils.h"
#include "dvec.h"
#include "waveletAMR.h"
namespace ot
{   


    /**
     * @brief Create a Mesh object by an array of octants. 
     * 
     * @param oct : pointer to a list of octant. 
     * @param num : number of 
     * @param eleOrder: element order. 
     * param comm: global communicator. 
     * @param verbose: if > 0 prints additional infomation on the Mesh generation. 
     * @param sm_type: scatter map type 
     * @param grain_sz : grain size for the partitoning, 
     * @param ld_tol : Load imbalance tolerance. 
     * @param sf_k : splitter fix values. 
     * @return Mesh mesh object. 
     */
    Mesh* createMesh(const ot::TreeNode* oct, unsigned int num,unsigned int eleOrder, MPI_Comm comm , unsigned int verbose=1, ot::SM_TYPE sm_type = ot::SM_TYPE::FDM, unsigned int grain_sz = DENDRO_DEFAULT_GRAIN_SZ, double ld_tol = DENDRO_DEFAULT_LB_TOL, unsigned int sf_k = DENDRO_DEFAULT_SF_K, unsigned int (*getWeight)(const ot::TreeNode *)=NULL);


    /**
     * @brief Generates an adaptive mesh based on the Wavelet AMR
     * 
     * @param func : spatially dependent function. 
     * @param wtol : wavelet tolerance. 
     * @param numVars : number of vars returned by the func. 
     * @param eleOrder : element order
     * @param comm: global communicator. 
     * @param refId : refinement var ids. 
     * @param verbose: if > 0 prints additional infomation on the Mesh generation. 
     * @param sm_type: scatter map type (FDM, FEM_CG, FEM_DG)
     * @param sz : number of refinement variables. 
     * @param grain_sz : grains size, 
     * @param ld_tol : load imbalance tolerance for flexible partitioning. 
     * @param sf_k : splitter fix k value (for large case runs)
     * @return Mesh* 
     */
    Mesh* createWAMRMesh(std::function<void(double,double,double,double*)> func, double wtol, unsigned int numVars, unsigned int eleOrder, MPI_Comm comm, unsigned int verbose=1,ot::SM_TYPE sm_type = ot::SM_TYPE::FDM, unsigned int * refIds=NULL, unsigned int sz=0, unsigned int grain_sz = DENDRO_DEFAULT_GRAIN_SZ, double ld_tol = DENDRO_DEFAULT_LB_TOL, unsigned int sf_k = DENDRO_DEFAULT_SF_K);

    /**
     * @brief Creates a mesh which is guranteed to converge 
     * @param pMesh : Pointer to current mesh object
     * @param wtol : wavelet tolerance, 
     * @param numVars : number of variables defined on the mesh. 
     * @param eleOrder : element order. 
     * @param refIds : refinment ids
     * @param sz : size of the refinement ids
     * @param maxiter : number of maximum iteration. 
     */
    void meshWAMRConvergence(ot::Mesh*& pMesh, std::function<void(double,double,double,double*)> func, double wtol, unsigned int numVars, unsigned int eleOrder,unsigned int * refIds=NULL, unsigned int sz=0, unsigned int maxiter=10);

    /**
     * @brief computes block unzip ghost node dependancies. 
     * 
     * @param pMesh : pointer to the mesh object
     * @param blk : blk id
     * @param gid : vector of ghost node ids related to the unzip, if the block is independent then the gid will have size zero. 
     */
    void computeBlockUnzipGhostNodes(const ot::Mesh* pMesh, unsigned int blk, std::vector<unsigned int>& gid);

    /**
     * @brief computes the block unzip padding elements. 
     * @param pMesh : pointer to the mesh object
     * @param blk : blk id
     * @param eid : vector of element ids. 
     */
    void computeBlockUnzipDepElements(const ot::Mesh* pMesh, unsigned int blk, std::vector<unsigned int>& eid);

    /**
     * @brief compute the corresponding time step level for a given block id. 
     * 
     * @param pMesh : pointer to mesh object. 
     * @param blk : local block id. 
     * @return unsigned int : corresponding time 
     */
    unsigned int  computeTLevel(const ot::Mesh* const pMesh, unsigned int blk);

    /**
     * @brief slice the mesh 
     * @param pMesh : mesh 
     * @param s_val : point on the slice plain
     * @param s_normal : normal vector to the slice. 
     * @param sids : elements on the slice
     * @return int 
     */
    int slice_mesh(const ot::Mesh* pMesh, unsigned int s_val[3], unsigned int s_normal[3], std::vector<unsigned int>& sids);

    /**
     * @brief Create a Split Mesh object with x left will have lmax and x right will have lmin level. 
     * 
     * @param eleOrder element order
     * @param lmin : level min
     * @param lmax : level max
     * @param comm : communicator mpi
     * @return Mesh* 
     */
    Mesh* createSplitMesh(unsigned int eleOrder, unsigned int lmin, unsigned int lmax,MPI_Comm comm);



}// end of namespae ot. 