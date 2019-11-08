/**
 * @file mathMeshUtils.h
 * @author: Milinda Fernando
 * @brief : Math utils using the mesh class. 
 * @version 0.1
 * @date 2019-11-07
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2019
 * 
 */

#pragma once

#include "mesh.h"

/**
 * @brief : Computes the volumetric norm $||v2-v1||_p$ with Reimann Sum.  
 * @note assumes the vectors have exchange the ghosts. 
 * @tparam T : type of the the vector
 * @param pMesh : Underlying mesh data 
 * @param v1 : input vector 1
 * @param v2 : intpur vector 2
 * @param p : order of the norm. 
 */
template<typename T>
double rsNormLp(const ot::Mesh* pMesh, const T* v1, const T* v2, unsigned int p)
{

    if(!pMesh->isActive())
        return -1.0;

    const unsigned int nPe = pMesh->getNumNodesPerElement();
    double * eleVec1 =new double[nPe];
    double * eleVec2 =new double[nPe];

    double rs=0;
    double rs_g=0;
    double hx;
    const ot::TreeNode* pNodes = pMesh->getAllElements().data();

    for(unsigned int e = pMesh->getElementLocalBegin(); e < pMesh->getElementLocalEnd(); e ++)
    {

        pMesh->getElementNodalValues(v1,eleVec1,e);
        pMesh->getElementNodalValues(v2,eleVec2,e);

        hx = ((1u<<(m_uiMaxDepth-pNodes[e].getLevel())))/(pMesh->getElementOrder());
        

        for(unsigned int n=0; n< nPe ; n++)
            rs += pow(fabs(eleVec1[n]-eleVec2[n]),p)*hx*hx*hx;
    }

    delete [] eleVec1;
    delete [] eleVec2;

    par::Mpi_Reduce(&rs,&rs_g,1,MPI_SUM,0,pMesh->getMPICommunicator());
    rs_g/= (1u<<(3*m_uiMaxDepth));
    rs_g = pow(rs_g,(1.0/p));
    return rs_g;

}