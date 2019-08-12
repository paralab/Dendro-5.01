//
// Created by milinda on 12/20/18.
//

/**
 * @brief contains ODA utility functions,
 * @author Milinda Fernando, School of Computing, University of Utah
 * */

#include "odaUtils.h"

void ot::computeODAFlags(const ot::Mesh *pMesh, std::vector<unsigned int> &flagList)
{
    flagList.clear();
    flagList.resize(pMesh->getAllElements().size(),0);

    if(!(pMesh->isActive())) return;

    const unsigned int * e2e=&(*(pMesh->getE2EMapping().begin()));
    const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));
    const unsigned int * e2n_dg=&(*(pMesh->getE2NMapping_DG().begin()));

    const unsigned int octantLocalBegin=pMesh->getElementLocalBegin();
    const unsigned int octantLocalEnd=pMesh->getElementLocalEnd();

    const unsigned int octantPreGhostBegin=pMesh->getElementPreGhostBegin();
    const unsigned int octantPreGhostEnd=pMesh->getElementPreGhostEnd();

    const unsigned int octantPostGhostBegin=pMesh->getElementPostGhostBegin();
    const unsigned int octantPostGhostEnd=pMesh->getElementPostGhostEnd();

    const unsigned int nodeLocalBegin=pMesh->getNodeLocalBegin();
    const unsigned int nodeLocalEnd=pMesh->getNodeLocalEnd();

    const unsigned int nodePreGhostBegin=pMesh->getNodePreGhostBegin();
    const unsigned int nodePreGhostEnd=pMesh->getNodePreGhostEnd();

    const unsigned int nodePostGhostBegin=pMesh->getNodePostGhostBegin();
    const unsigned int nodePostGhostEnd=pMesh->getNodePostGhostEnd();

    const unsigned int numDirections=pMesh->getNumDirections();

    const unsigned int nPe=pMesh->getNumNodesPerElement();
    const ot::TreeNode* allElements=&(*(pMesh->getAllElements().begin()));

    const unsigned int octTreeMin=0;
    const unsigned int octTreeMax=(1u<<m_uiMaxDepth);

    unsigned int lookUp;
    bool isWDep=false;


    // 1. flag the  local elements.
    for(unsigned int ele=octantLocalBegin;ele<octantLocalEnd;ele++)
    {
        isWDep= false;
        for(unsigned int node=0;node<nPe;node++)
        {
            lookUp=e2n[ele*nPe+node];
            if((lookUp!=LOOK_UP_TABLE_DEFAULT) && !(lookUp>=nodeLocalBegin && lookUp<nodeLocalEnd))
            {
                binOp::setBit(flagList[ele],ODA_W_DEPENDENT_FLAG_BIT);
                isWDep=true;
                break;
            }
        }

        if(!isWDep)
        {
            binOp::setBit(flagList[ele],ODA_INDEPENDENT_FLAG_BIT);

        }



        if((allElements[ele].minX()==octTreeMin || allElements[ele].minY()==octTreeMin || allElements[ele].minZ()==octTreeMin
            || allElements[ele].maxX()==octTreeMax || allElements[ele].maxY()==octTreeMax || allElements[ele].maxZ() == octTreeMax) )
        {
            binOp::setBit(flagList[ele],ODA_W_BOUNDARY_FLAG_BIT);
        }

    }
    // do not use this this is only the ghost element exchanged at round 1. 
    //const std::vector<unsigned int >& lev1_Octants=pMesh->getLevel1GhostElementIndices();

    for(unsigned int ele=octantPreGhostBegin;ele<octantPreGhostEnd;ele++)
    {

        //binOp::setBit(m_uiOctantFlags[lev1_Octants[g]],ODA_W_DEPENDENT_FLAG_BIT);
        for(unsigned int node=0;node<nPe;node++)
        {
            lookUp=e2n[ele*nPe+node];
            if((lookUp>=nodeLocalBegin && lookUp<nodeLocalEnd))
            {
                binOp::setBit(flagList[ele],ODA_W_DEPENDENT_FLAG_BIT);
                break;
            }

        }
    }


    for(unsigned int ele=octantPostGhostBegin;ele<octantPostGhostEnd;ele++)
    {

        //binOp::setBit(m_uiOctantFlags[lev1_Octants[g]],ODA_W_DEPENDENT_FLAG_BIT);
        for(unsigned int node=0;node<nPe;node++)
        {
            lookUp=e2n[ele*nPe+node];
            if((lookUp>=nodeLocalBegin && lookUp<nodeLocalEnd))
            {
                binOp::setBit(flagList[ele],ODA_W_DEPENDENT_FLAG_BIT);
                break;
            }

        }
    }




}


void ot::computeLocalToGlobalNodalMap(const ot::Mesh* pMesh,std::vector<DendroIntL>& map,DendroIntL& globalNodeSz)
{

    map.clear();
    
    if(pMesh->isActive())
    {
        const unsigned int npesAll=pMesh->getMPICommSizeGlobal();
        pMesh->createVector(map);

        const unsigned int nodeLocalBegin=pMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd=pMesh->getNodeLocalEnd();

        const unsigned int nodePreGhostBegin=pMesh->getNodePreGhostBegin();
        const unsigned int nodePreGhostEnd=pMesh->getNodePreGhostEnd();

        const unsigned int nodePostGhostBegin=pMesh->getNodePostGhostBegin();
        const unsigned int nodePostGhostEnd=pMesh->getNodePostGhostEnd();

        if(npesAll==1) {

            for(unsigned int node=nodeLocalBegin;node<nodeLocalEnd;node++)
                map[node]=node;
            
            globalNodeSz=pMesh->getNumLocalMeshNodes();

        }else{

            for(unsigned int node=nodePreGhostBegin;node<nodePreGhostEnd;node++)
                map[node]=0;

            for(unsigned int node=nodePostGhostBegin;node<nodePostGhostEnd;node++)
                map[node]=0;

            const unsigned int activeNpes=pMesh->getMPICommSize();
            const unsigned int activeRank=pMesh->getMPIRank();

            MPI_Comm activeComm=pMesh->getMPICommunicator();

            unsigned int localNodes=pMesh->getNumLocalMeshNodes();
            unsigned int * nodalCount=new unsigned int [activeNpes];
            unsigned int * nodalOffset=new unsigned int [activeNpes];

            for(unsigned int p=0;p<activeNpes;p++)
                nodalCount[p]=0;

            par::Mpi_Allgather(&localNodes,nodalCount,1,activeComm);

            nodalOffset[0]=0;
            omp_par::scan(nodalCount,nodalOffset,activeNpes);

            globalNodeSz=nodalOffset[activeNpes-1] + nodalCount[activeNpes-1];

            for(unsigned int node=nodeLocalBegin;node<nodeLocalEnd;node++)
                map[node]=nodalOffset[activeRank] + (node-nodeLocalBegin) ;


            delete [] nodalCount;
            delete [] nodalOffset;

            const unsigned int * sendNodeOffset = &(*(pMesh->getNodalSendOffsets().begin()));
            const unsigned int * sendNodeCounts = &(*(pMesh->getNodalSendCounts().begin()));

            const unsigned int * recvNodeOffset = &(*(pMesh->getNodalRecvOffsets().begin()));
            const unsigned int * recvNodeCounts = &(*(pMesh->getNodalRecvCounts().begin()));

            const unsigned int * sendNodeSM = &(*(pMesh->getSendNodeSM().begin()));
            const unsigned int * recvNodeSM = &(*(pMesh->getRecvNodeSM().begin()));


            DendroIntL* sendBufferNodes = new DendroIntL[sendNodeOffset[activeNpes-1] + sendNodeCounts[activeNpes-1]];
            DendroIntL* recvBufferNodes = new DendroIntL[recvNodeOffset[activeNpes-1] + recvNodeCounts[activeNpes-1]];


            for(unsigned int p=0;p<activeNpes;p++) {
                for (unsigned int k = sendNodeOffset[p]; k < (sendNodeOffset[p] + sendNodeCounts[p]); k++) {
                    sendBufferNodes[k] = map[sendNodeSM[k]];
                }
            }

    #ifdef ALLTOALL_SPARSE
            par::Mpi_Alltoallv_sparse(sendBufferNodes,(int *)sendNodeCounts,(int *)sendNodeOffset,recvBufferNodes,(int *) recvNodeCounts,(int *)recvNodeOffset,activeComm);
    #else
            par::Mpi_Alltoallv(sendBufferNodes,(int *)sendNodeCounts,(int *)sendNodeOffset,recvBufferNodes,(int *) recvNodeCounts,(int *)recvNodeOffset,activeComm);
    #endif

            for(unsigned int p=0;p<activeNpes;p++)
            {
                for(unsigned int k=recvNodeOffset[p];k<(recvNodeOffset[p]+recvNodeCounts[p]);k++)
                {
                    map[recvNodeSM[k]]=recvBufferNodes[k];
                }
            }


            delete [] sendBufferNodes;
            delete [] recvBufferNodes;
        }
        
    }
    
    par::Mpi_Bcast(&globalNodeSz,1,0,pMesh->getMPIGlobalCommunicator());
    
    









}
