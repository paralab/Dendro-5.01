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
#include <functional>
#include <iostream>
#include <vector>

#include "TreeNode.h"
#include "asyncExchangeContex.h"
#include "dvec.h"
#include "mesh.h"
#include "octUtils.h"
#include "parUtils.h"
#include "waveletAMR.h"
namespace ot {

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
Mesh* createMesh(const ot::TreeNode* oct, unsigned int num,
                 unsigned int eleOrder, MPI_Comm comm, unsigned int verbose = 1,
                 ot::SM_TYPE sm_type   = ot::SM_TYPE::FDM,
                 unsigned int grain_sz = DENDRO_DEFAULT_GRAIN_SZ,
                 double ld_tol         = DENDRO_DEFAULT_LB_TOL,
                 unsigned int sf_k     = DENDRO_DEFAULT_SF_K,
                 unsigned int (*getWeight)(const ot::TreeNode*) = NULL);

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
Mesh* createWAMRMesh(std::function<void(double, double, double, double*)> func,
                     double wtol, unsigned int numVars, unsigned int eleOrder,
                     MPI_Comm comm, unsigned int verbose = 1,
                     ot::SM_TYPE sm_type  = ot::SM_TYPE::FDM,
                     unsigned int* refIds = NULL, unsigned int sz = 0,
                     unsigned int grain_sz = DENDRO_DEFAULT_GRAIN_SZ,
                     double ld_tol         = DENDRO_DEFAULT_LB_TOL,
                     unsigned int sf_k     = DENDRO_DEFAULT_SF_K);

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
void meshWAMRConvergence(
    ot::Mesh*& pMesh, std::function<void(double, double, double, double*)> func,
    double wtol, unsigned int numVars, unsigned int eleOrder,
    unsigned int* refIds = NULL, unsigned int sz = 0,
    unsigned int maxiter = 10);

/**
 * @brief computes block unzip ghost node dependancies.
 *
 * @param pMesh : pointer to the mesh object
 * @param blk : blk id
 * @param gid : vector of ghost node ids related to the unzip, if the block is
 * independent then the gid will have size zero.
 */
void computeBlockUnzipGhostNodes(const ot::Mesh* pMesh, unsigned int blk,
                                 std::vector<unsigned int>& gid);

/**
 * @brief computes the block unzip padding elements.
 * @param pMesh : pointer to the mesh object
 * @param blk : blk id
 * @param eid : vector of element ids.
 */
void computeBlockUnzipDepElements(const ot::Mesh* pMesh, unsigned int blk,
                                  std::vector<unsigned int>& eid);

/**
 * @brief compute the corresponding time step level for a given block id.
 *
 * @param pMesh : pointer to mesh object.
 * @param blk : local block id.
 * @return unsigned int : corresponding time
 */
unsigned int computeTLevel(const ot::Mesh* const pMesh, unsigned int blk);

/**
 * @brief slice the mesh
 * @param pMesh : mesh
 * @param s_val : point on the slice plain
 * @param s_normal : normal vector to the slice.
 * @param sids : elements on the slice
 * @return int
 */
int slice_mesh(const ot::Mesh* pMesh, unsigned int s_val[3],
               unsigned int s_normal[3], std::vector<unsigned int>& sids);

/**
 * @brief Create a Split Mesh object with x left will have lmax and x right will
 * have lmin level.
 *
 * @param eleOrder element order
 * @param lmin : level min
 * @param lmax : level max
 * @param comm : communicator mpi
 * @return Mesh*
 */
Mesh* createSplitMesh(unsigned int eleOrder, unsigned int lmin,
                      unsigned int lmax, MPI_Comm comm);

template <typename T>
double computeElementIntegral(T* elementalVec, const ot::TreeNode& node,
                              const RefElement* refEl, T& elementIntegral,
                              const T Lx, const T Ly, const T Lz) {
    // this will compute the elemental integration
    const unsigned int N         = refEl->get1DNumInterpolationPoints();
    // interpolation *and* weighting matrix - calculated W @ Q1d
    const double* w_Q1d          = refEl->getWMulQ1D();

    const unsigned int len       = (1u << (m_uiMaxDepth - node.getLevel()));
    // dx, dy, and dz are all then just the length divided by the number of
    // points -1, and are the same as h
    // const double h         = (double)len / (double)refEl->getElementOrder();
    const double full_oct_domain = (double)(1u << m_uiMaxDepth);

    // calculate the integral with the scaled/transformed weight vector
    elementIntegral              = 0.0;
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                elementIntegral += w_Q1d[i] * w_Q1d[j] * w_Q1d[k] *
                                   elementalVec[k * N * N + j * N + i];
            }
        }
    }

    // with h now calculated, we need to scale by the actual block length,
    // i.e. jacobian = (phys_x_len) * (phys_y_len) * (phys_z_len) / 8
    const double Lx_scale = len * (Lx / full_oct_domain);
    const double Ly_scale = len * (Ly / full_oct_domain);
    const double Lz_scale = len * (Lz / full_oct_domain);
    const double jacobian = (Lx_scale * Ly_scale * Lz_scale) / 8.0;
    elementIntegral *= jacobian;

    return elementIntegral;
}

template <typename T>
T integrateScalarFieldOnMesh(ot::Mesh* mesh, T* in) {
    // just return 0.0 if inactive mesh
    if (!mesh->isActive()) return 0.0;

    T localIntegral = 0.0;
    // perform ghost exchange to sync values
    mesh->performGhostExchange(in);

    const unsigned int nPe      = mesh->getNumNodesPerElement();
    const unsigned int eleOrder = mesh->getElementOrder();

    T* elementalVec             = new T[nPe];
    T elementIntegral;

    const T Lx = mesh->getDomainMaxPt().x() - mesh->getDomainMinPt().x();
    const T Ly = mesh->getDomainMaxPt().y() - mesh->getDomainMinPt().y();
    const T Lz = mesh->getDomainMaxPt().z() - mesh->getDomainMinPt().z();

    // process the ghost elements for shared faces and edges
    for (unsigned int ele = mesh->getElementLocalBegin();
         ele < mesh->getElementLocalEnd(); ele++) {
        mesh->getElementNodalValues(in, elementalVec, ele);

        // then compute the element-based integral
        computeElementIntegral(elementalVec, mesh->getAllElements()[ele],
                               mesh->getReferenceElement(), elementIntegral, Lx,
                               Ly, Lz);

        localIntegral += elementIntegral;
    }

    delete[] elementalVec;

    // then all reduce to root
    T globalIntegral = 0.0;
    par::Mpi_Allreduce(&localIntegral, &globalIntegral, 1, MPI_SUM,
                       mesh->getMPICommunicator());

    return globalIntegral;
}

template <typename T>
T calculateL2ErrorFullMeshIntegration(ot::Mesh* mesh, T* in, T* truth) {
    if (!mesh->isActive()) return 0.0;

    // create a new vector that's the same size as the in and truth, since we're
    // using mesh we'll have the same size as a "zipped" vector
    std::vector<T> squaredDiffVector;
    mesh->createVector(squaredDiffVector);

    // then iterate through the values
    for (size_t i = 0; i < mesh->getDegOfFreedom(); i++) {
        squaredDiffVector[i] = (in[i] - truth[i]) * (in[i] - truth[i]);
    }

    // NOTE: performGhostExchange happens inside this function
    return integrateScalarFieldOnMesh(mesh, squaredDiffVector.data());
}

template <typename T>
T calculateL2FullMeshIntegration(ot::Mesh* mesh, T* in) {
    if (!mesh->isActive()) return 0.0;

    // create a new vector that's the same size as the in, since we're
    // using mesh we'll have the same size as a "zipped" vector
    std::vector<T> squaredVector;
    mesh->createVector(squaredVector);

    // then iterate through the values
    for (size_t i = 0; i < mesh->getDegOfFreedom(); i++) {
        squaredVector[i] = (in[i]) * (in[i]);
    }

    // NOTE: performGhostExchange happens inside this function
    return integrateScalarFieldOnMesh(mesh, squaredVector.data());
}

template <typename T>
void alloc_mpi_ctx(const Mesh* pMesh,
                   std::vector<AsyncExchangeContex>& ctx_list, int dof,
                   int async_k) {
    if (pMesh->getMPICommSizeGlobal() == 1 || !pMesh->isActive()) return;

    {
        const std::vector<unsigned int>& nodeSendCount =
            pMesh->getNodalSendCounts();
        const std::vector<unsigned int>& nodeSendOffset =
            pMesh->getNodalSendOffsets();

        const std::vector<unsigned int>& e_sf = pMesh->getElementSendOffsets();
        const std::vector<unsigned int>& e_sc = pMesh->getElementSendCounts();
        const std::vector<unsigned int>& e_rf = pMesh->getElementRecvOffsets();
        const std::vector<unsigned int>& e_rc = pMesh->getElementRecvCounts();

        const std::vector<unsigned int>& nodeRecvCount =
            pMesh->getNodalRecvCounts();
        const std::vector<unsigned int>& nodeRecvOffset =
            pMesh->getNodalRecvOffsets();

        const std::vector<unsigned int>& sendProcList =
            pMesh->getSendProcList();
        const std::vector<unsigned int>& recvProcList =
            pMesh->getRecvProcList();

        const std::vector<unsigned int>& sendNodeSM = pMesh->getSendNodeSM();
        const std::vector<unsigned int>& recvNodeSM = pMesh->getRecvNodeSM();

        const unsigned int activeNpes               = pMesh->getMPICommSize();
        const unsigned int nPe = pMesh->getNumNodesPerElement();

        const unsigned int sendBSzCg =
            (nodeSendOffset[activeNpes - 1] + nodeSendCount[activeNpes - 1]);
        const unsigned int recvBSzCg =
            (nodeRecvOffset[activeNpes - 1] + nodeRecvCount[activeNpes - 1]);

        const unsigned int sendBSzDg =
            (e_sf[activeNpes - 1] + e_sc[activeNpes - 1]) * nPe;
        const unsigned int recvBSzDg =
            (e_rf[activeNpes - 1] + e_rc[activeNpes - 1]) * nPe;

        const unsigned int sendBSz = std::max(sendBSzCg, sendBSzDg);
        const unsigned int recvBSz = std::max(recvBSzCg, recvBSzDg);

        ctx_list.resize(async_k);
        for (unsigned int i = 0; i < async_k; i++) {
            const unsigned int v_begin  = ((i * dof) / async_k);
            const unsigned int v_end    = (((i + 1) * dof) / async_k);
            const unsigned int batch_sz = (v_end - v_begin);

            if (sendBSz)
                ctx_list[i].allocateSendBuffer(batch_sz * sendBSz * sizeof(T));

            if (recvBSz)
                ctx_list[i].allocateRecvBuffer(batch_sz * recvBSz * sizeof(T));

            ctx_list[i].m_send_req.resize(pMesh->getMPICommSize(),
                                          MPI_Request());
            ctx_list[i].m_recv_req.resize(pMesh->getMPICommSize(),
                                          MPI_Request());
        }
    }

    return;
}

template <typename T>
void dealloc_mpi_ctx(const Mesh* pMesh,
                     std::vector<AsyncExchangeContex>& ctx_list, int dof,
                     int async_k) {
    if (pMesh->getMPICommSizeGlobal() == 1 || !pMesh->isActive()) return;

    const std::vector<unsigned int>& nodeSendCount =
        pMesh->getNodalSendCounts();
    const std::vector<unsigned int>& nodeSendOffset =
        pMesh->getNodalSendOffsets();

    const std::vector<unsigned int>& nodeRecvCount =
        pMesh->getNodalRecvCounts();
    const std::vector<unsigned int>& nodeRecvOffset =
        pMesh->getNodalRecvOffsets();

    const std::vector<unsigned int>& sendProcList = pMesh->getSendProcList();
    const std::vector<unsigned int>& recvProcList = pMesh->getRecvProcList();

    const std::vector<unsigned int>& sendNodeSM   = pMesh->getSendNodeSM();
    const std::vector<unsigned int>& recvNodeSM   = pMesh->getRecvNodeSM();

    const unsigned int activeNpes                 = pMesh->getMPICommSize();

    const unsigned int sendBSz =
        nodeSendOffset[activeNpes - 1] + nodeSendCount[activeNpes - 1];
    const unsigned int recvBSz =
        nodeRecvOffset[activeNpes - 1] + nodeRecvCount[activeNpes - 1];

    for (unsigned int i = 0; i < ctx_list.size(); i++) {
        const unsigned int v_begin  = ((i * dof) / async_k);
        const unsigned int v_end    = (((i + 1) * dof) / async_k);
        const unsigned int batch_sz = (v_end - v_begin);

        if (sendBSz) ctx_list[i].deAllocateSendBuffer();

        if (recvBSz) ctx_list[i].deAllocateRecvBuffer();

        ctx_list[i].m_send_req.clear();
        ctx_list[i].m_recv_req.clear();
    }

    ctx_list.clear();
}

}  // namespace ot
