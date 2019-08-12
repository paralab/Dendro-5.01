/**
 * @file nuts.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Spatially adaptive non uniform time stepping framework. 
 * @version 0.1
 * @date 2019-07-12
 * @copyright Copyright (c) 2019, School of Computing, University of Utah. 
 * 
 */

#include "nuts.h"

namespace ts
{
    NUTS::NUTS(ot::Mesh* octMesh, unsigned int s, double tb, double te, Point ptMin, Point ptMax, double cfl)
    {
        m_uiOctGrid = octMesh;
        m_uiS = s;
        m_uiTB = tb;
        m_uiTE = te;
        m_uiCFL = cfl;
        isActive = m_uiOctGrid->isActive();
        m_uiCommGlobal = m_uiOctGrid->getMPIGlobalCommunicator();

        m_uiPtMin = ptMin;
        m_uiPtMax = ptMax;

        const double dx = ( ptMax.x() - ptMin.x() ) / (1u<<m_uiMaxDepth);
        const double dy = ( ptMax.y() - ptMin.y() ) / (1u<<m_uiMaxDepth);
        const double dz = ( ptMax.z() - ptMin.z() ) / (1u<<m_uiMaxDepth);
        
        if(isActive)
            m_uiCommActive = m_uiOctGrid->getMPICommunicator();
        else
            m_uiCommActive = MPI_COMM_NULL;

        unsigned int activeNpes;
        std::vector<ot::TreeNode> blkNode;

        // identify dependent and independent blocks. 
        m_uiOctGrid->flagBlocks();

        if(isActive)
        {

            const ot::TreeNode * allElements = &(*(m_uiOctGrid->getAllElements().begin()));
            const std::vector<ot::Block> blkList = m_uiOctGrid->getLocalBlockList();

            
            blkNode.resize(blkList.size());
            for(unsigned int i=0;i<blkList.size();i++)
                blkNode[i] = blkList[i].getBlockNode();
            

            unsigned int l_min,l_max,l_min_g, l_max_g;
            const unsigned int numAllElements = m_uiOctGrid->getAllElements().size();
            const unsigned int localEleBegin = m_uiOctGrid->getElementLocalBegin();
            const unsigned int localEleEnd = m_uiOctGrid->getElementLocalEnd();

            l_min = allElements[localEleBegin].getLevel();
            l_max = allElements[localEleBegin].getLevel();

            for(unsigned int ele = localEleBegin; ele< localEleEnd; ele++)
            {
                if(l_min > allElements[ele].getLevel())
                    l_min = allElements[ele].getLevel();
                
                if(l_max < allElements[ele].getLevel())
                    l_max = allElements[ele].getLevel();

            }

            par::Mpi_Allreduce(&l_min,&l_min_g,1,MPI_MIN,m_uiCommActive);
            par::Mpi_Allreduce(&l_max,&l_max_g,1,MPI_MAX,m_uiCommActive);

            // compute the dt min and max. 

            m_uiDTMin = m_uiCFL * (1u<<l_min) * dx;
            m_uiDTMax = m_uiCFL * (1u<<l_max) * dx;
            
            activeNpes = m_uiOctGrid->getMPICommSize();

        }
        
        par::Mpi_Bcast(&m_uiDTMin,1,0,m_uiCommGlobal);
        par::Mpi_Bcast(&m_uiDTMax,1,0,m_uiCommGlobal);
        par::Mpi_Bcast(&activeNpes,1,0,m_uiCommGlobal);
        
        m_uiBlockGrid = new ot::Mesh(blkNode,1,1,m_uiCommGlobal,false,ot::FDM);

        
    }

    NUTS::~NUTS()
    {
        delete m_uiBlockGrid;
    }

    void NUTS::computeBlockSM()
    {

        if(m_uiOctGrid->isActive())
        {
            const ot::TreeNode * pNodes= &(*(m_uiOctGrid->getAllElements().begin()));

            const unsigned int blockLocalBegin = m_uiBlockGrid->getElementLocalBegin();
            const unsigned int blockLocalEnd  = m_uiBlockGrid-> getElementLocalEnd();

            const unsigned int nodeLocalBegin = m_uiOctGrid->getNodeLocalBegin();
            const unsigned int nodeLocalEnd = m_uiOctGrid->getNodeLocalEnd();

            const unsigned int* e2n_cg = &(*(m_uiOctGrid->getE2NMapping().begin()));
            const unsigned int* e2e = &(*(m_uiOctGrid->getE2EMapping().begin()));

            const std::vector<ot::Block> blkList = m_uiOctGrid->getLocalBlockList();
            const unsigned int nPe = m_uiOctGrid->getNumNodesPerElement();

            unsigned int lookup,node_cg;
            std::vector<unsigned int> * ghostNodesByLevel  = new std::vector<unsigned int>[m_uiMaxDepth];
            unsigned int child[NUM_CHILDREN];

            for(unsigned int blk=0; blk< blkList.size(); blk++)
            {
                if(blkList[blk].getBlockType() != ot::BlockType::FULLY_INDEPENDENT)
                {   
                    // these block depends on some ghost nodal values. 
                    const unsigned int pWidth = blkList[blk].get1DPadWidth();
                    const ot::TreeNode blkNode = blkList[blk].getBlockNode();
                    const unsigned int regLevel = blkList[blk].getRegularGridLev();

                    for(unsigned int elem = blkList[blk].getLocalElementBegin(); elem < blkList[blk].getLocalElementEnd(); elem++)
                    {
                        const unsigned int ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLevel);
                        const unsigned int ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLevel);
                        const unsigned int ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLevel);

                        const unsigned int emin = 0;
                        const unsigned int emax = (1u<<(regLevel-blkNode.getLevel()))-1;

                        if(((ei == emin) || (ej == emin) || (ek== emin) || (ei==emax) || (ej==emax) || (ek==emax)))
                        {
                            // this is block boundary element. 
                            for(unsigned int node =0; node < nPe; node++)
                            {
                                node_cg = e2n_cg[elem*nPe + node];
                                if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                    ghostNodesByLevel[regLevel].push_back(node_cg); 
                            }


                        }

                        if(pWidth > 0)
                        {   
                            // we need to look for the boundary neigbours only when the padding width is > 0 . 
                            if(ei==emin)
                            { 
                                // OCT_DIR_LEFT
                                lookup = e2e[elem*NUM_FACES + OCT_DIR_LEFT];
                                if(lookup!=LOOK_UP_TABLE_DEFAULT)
                                {
                                    if(pNodes[lookup].getLevel() > regLevel )
                                    {  
                                       // neighbour octant is smaller. 
                                       assert(pNodes[lookup].getLevel()==(regLevel+1));
                                       //child.resize(NUM_CHILDREN,LOOK_UP_TABLE_DEFAULT);
                                       // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                                       child[1]=lookup;
                                       child[3]=e2e[child[1]*NUM_FACES+OCT_DIR_UP];
                                       assert(child[3]!=LOOK_UP_TABLE_DEFAULT);
                                       child[5]=e2e[child[1]*NUM_FACES+OCT_DIR_FRONT];
                                       assert(child[5]!=LOOK_UP_TABLE_DEFAULT);
                                       child[7]=e2e[child[3]*NUM_FACES+OCT_DIR_FRONT];
                                       assert(child[7]!=LOOK_UP_TABLE_DEFAULT);

                                       child[0]=LOOK_UP_TABLE_DEFAULT; 
                                       child[2]=LOOK_UP_TABLE_DEFAULT; 
                                       child[4]=LOOK_UP_TABLE_DEFAULT; 
                                       child[6]=LOOK_UP_TABLE_DEFAULT; 

                                       for(unsigned int c=0;c<NUM_CHILDREN; c++)
                                       {
                                           if( child[c] != LOOK_UP_TABLE_DEFAULT )
                                           {
                                               for(unsigned int node =0; node < nPe; node++)
                                               {
                                                   node_cg = e2n_cg[child[c]*nPe + node];
                                                   if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                    ghostNodesByLevel[regLevel].push_back(node_cg);
                                               }
                                       
                                            }
                                       }
                                
                                    }else
                                    {   
                                        // neighbour octant is same lev or coarser
                                        assert(pNodes[lookup].getLevel() <= regLevel );
                                        for(unsigned int node = 0; node < nPe; node++)
                                        {
                                            node_cg = e2n_cg[lookup * nPe + node];
                                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                ghostNodesByLevel[regLevel].push_back(node_cg);
                                        }


                                    }

                                }

                                    
                            }

                            if(ei==emax)
                            { // OCT_DIR_RIGHT

                                lookup = e2e[elem*NUM_FACES + OCT_DIR_RIGHT];
                                if(lookup!=LOOK_UP_TABLE_DEFAULT)
                                {
                                    if(pNodes[lookup].getLevel() > regLevel )
                                    { // neighbour octant is smaller. 
                                      
                                      // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                                      child[0]=lookup;
                                      child[2]=e2e[child[0]*NUM_FACES+OCT_DIR_UP];
                                      assert(child[2]!=LOOK_UP_TABLE_DEFAULT);
                                      child[4]=e2e[child[0]*NUM_FACES+OCT_DIR_FRONT];
                                      assert(child[4]!=LOOK_UP_TABLE_DEFAULT);
                                      child[6]=e2e[child[2]*NUM_FACES+OCT_DIR_FRONT];
                                      assert(child[6]!=LOOK_UP_TABLE_DEFAULT);

                                      child[1]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_RIGHT];
                                      child[3]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_RIGHT];
                                      child[5]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_RIGHT];
                                      child[7]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_RIGHT];

                                      for(unsigned int c=0;c<NUM_CHILDREN; c++)
                                       {
                                           if( child[c] != LOOK_UP_TABLE_DEFAULT )
                                           {
                                               for(unsigned int node =0; node < nPe; node++)
                                               {
                                                   node_cg = e2n_cg[child[c]*nPe + node];
                                                   if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                    ghostNodesByLevel[regLevel].push_back(node_cg);
                                               }
                                       
                                            }
                                       }
                                      

                                    }else
                                    {   
                                        // neighbour octant is same lev or coarser
                                        assert(pNodes[lookup].getLevel() <= regLevel );
                                        for(unsigned int node = 0; node < nPe; node++)
                                        {
                                            node_cg = e2n_cg[lookup * nPe + node];
                                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                ghostNodesByLevel[regLevel].push_back(node_cg);
                                        }

                                    }

                                }
                            }

                            if(ej==emin)
                            {   // OCT_DIR_DOWN
                                lookup = e2e[elem*NUM_FACES + OCT_DIR_DOWN];
                                if(lookup!=LOOK_UP_TABLE_DEFAULT)
                                {

                                    if(pNodes[lookup].getLevel() > regLevel )
                                    { // neighbour octant is smaller. 

                                        child[2]=lookup;
                                        child[3]=e2e[child[2]*NUM_FACES+OCT_DIR_RIGHT];
                                        assert(child[3]!=LOOK_UP_TABLE_DEFAULT);
                                        child[6]=e2e[child[2]*NUM_FACES+OCT_DIR_FRONT];
                                        assert(child[6]!=LOOK_UP_TABLE_DEFAULT);
                                        child[7]=e2e[child[3]*NUM_FACES+OCT_DIR_FRONT];
                                        assert(child[7]!=LOOK_UP_TABLE_DEFAULT);

                                        child[0]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_DOWN];
                                        child[1]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_DOWN];
                                        child[4]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_DOWN];
                                        child[5]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_DOWN];

                                        for(unsigned int c=0;c<NUM_CHILDREN; c++)
                                        {
                                            if( child[c] != LOOK_UP_TABLE_DEFAULT )
                                            {
                                                for(unsigned int node =0; node < nPe; node++)
                                                {
                                                    node_cg = e2n_cg[child[c]*nPe + node];
                                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                        ghostNodesByLevel[regLevel].push_back(node_cg);
                                                }
                                        
                                            }
                                        }

                                    }else
                                    { // neighbour octant is same lev or coarser
                                        
                                        assert(pNodes[lookup].getLevel() <= regLevel );
                                        for(unsigned int node = 0; node < nPe; node++)
                                        {
                                            node_cg = e2n_cg[lookup * nPe + node];
                                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                ghostNodesByLevel[regLevel].push_back(node_cg);
                                        }

                                    }

                                }
                            }

                            if(ej==emax)
                            {   // OCT_DIR_UP
                                lookup = e2e[elem*NUM_FACES + OCT_DIR_UP];
                                if(lookup!=LOOK_UP_TABLE_DEFAULT)
                                {
                                   if(pNodes[lookup].getLevel() > regLevel )
                                    { // neighbour octant is smaller. 

                                        // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                                        child[0]=lookup;
                                        child[1]=e2e[child[0]*NUM_FACES+OCT_DIR_RIGHT];
                                        assert(child[1]!=LOOK_UP_TABLE_DEFAULT);
                                        child[4]=e2e[child[0]*NUM_FACES+OCT_DIR_FRONT];
                                        assert(child[4]!=LOOK_UP_TABLE_DEFAULT);
                                        child[5]=e2e[child[1]*NUM_FACES+OCT_DIR_FRONT];
                                        assert(child[5]!=LOOK_UP_TABLE_DEFAULT);

                                        child[2]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_UP];
                                        child[3]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_UP];
                                        child[6]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_UP];
                                        child[7]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_UP];

                                        for(unsigned int c=0;c<NUM_CHILDREN; c++)
                                        {
                                            if( child[c] != LOOK_UP_TABLE_DEFAULT )
                                            {
                                                for(unsigned int node =0; node < nPe; node++)
                                                {
                                                    node_cg = e2n_cg[child[c]*nPe + node];
                                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                        ghostNodesByLevel[regLevel].push_back(node_cg);
                                                }
                                        
                                                }
                                        }

                                    }else
                                    { // neighbour octant is same lev or coarser

                                        assert(pNodes[lookup].getLevel() <= regLevel );
                                        for(unsigned int node = 0; node < nPe; node++)
                                        {
                                            node_cg = e2n_cg[lookup * nPe + node];
                                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                ghostNodesByLevel[regLevel].push_back(node_cg);
                                        }

                                    }
                                }
                            }


                            if(ek==emin)
                            {   // OCT_DIR_BACK
                                lookup = e2e[elem*NUM_FACES + OCT_DIR_BACK];
                                if(lookup!=LOOK_UP_TABLE_DEFAULT)
                                {
                                    if(pNodes[lookup].getLevel() > regLevel )
                                    { // neighbour octant is smaller. 

                                        child[4]=lookup;
                                        child[5]=e2e[child[4]*NUM_FACES+OCT_DIR_RIGHT];
                                        assert(child[5]!=LOOK_UP_TABLE_DEFAULT);
                                        child[6]=e2e[child[4]*NUM_FACES+OCT_DIR_UP];
                                        assert(child[6]!=LOOK_UP_TABLE_DEFAULT);
                                        child[7]=e2e[child[5]*NUM_FACES+OCT_DIR_UP];
                                        assert(child[7]!=LOOK_UP_TABLE_DEFAULT);

                                        child[0]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_BACK];
                                        child[1]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_BACK];
                                        child[2]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_BACK];
                                        child[3]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_BACK];

                                        for(unsigned int c=0;c<NUM_CHILDREN; c++)
                                        {
                                            if( child[c] != LOOK_UP_TABLE_DEFAULT )
                                            {
                                                for(unsigned int node =0; node < nPe; node++)
                                                {
                                                    node_cg = e2n_cg[child[c]*nPe + node];
                                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                        ghostNodesByLevel[regLevel].push_back(node_cg);
                                                }
                                        
                                                }
                                        }

                                    }else
                                    { // neighbour octant is same lev or coarser

                                        assert(pNodes[lookup].getLevel() <= regLevel );
                                        for(unsigned int node = 0; node < nPe; node++)
                                        {
                                            node_cg = e2n_cg[lookup * nPe + node];
                                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                ghostNodesByLevel[regLevel].push_back(node_cg);
                                        }

                                    }


                                }
                            }

                            if(ek==emax)
                            {  // OCT_DIR_FRONT
                                lookup = e2e[elem*NUM_FACES + OCT_DIR_FRONT];
                                if(lookup!=LOOK_UP_TABLE_DEFAULT)
                                {
                                    if(pNodes[lookup].getLevel() > regLevel )
                                    { // neighbour octant is smaller. 

                                        child[0]=lookup;
                                        child[1]=e2e[child[0]*NUM_FACES+OCT_DIR_RIGHT];
                                        assert(child[1]!=LOOK_UP_TABLE_DEFAULT);
                                        child[2]=e2e[child[0]*NUM_FACES+OCT_DIR_UP];
                                        assert(child[2]!=LOOK_UP_TABLE_DEFAULT);
                                        child[3]=e2e[child[1]*NUM_FACES+OCT_DIR_UP];
                                        assert(child[3]!=LOOK_UP_TABLE_DEFAULT);

                                        child[4]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_FRONT];
                                        child[5]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_FRONT];
                                        child[6]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_FRONT];
                                        child[7]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_FRONT];

                                        for(unsigned int c=0;c<NUM_CHILDREN; c++)
                                        {
                                            if( child[c] != LOOK_UP_TABLE_DEFAULT )
                                            {
                                                for(unsigned int node =0; node < nPe; node++)
                                                {
                                                    node_cg = e2n_cg[child[c]*nPe + node];
                                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                        ghostNodesByLevel[regLevel].push_back(node_cg);
                                                }
                                        
                                            }
                                        }

                                    }else
                                    { // neighbour octant is same lev or coarser

                                        assert(pNodes[lookup].getLevel() <= regLevel );
                                        for(unsigned int node = 0; node < nPe; node++)
                                        {
                                            node_cg = e2n_cg[lookup * nPe + node];
                                            if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                                ghostNodesByLevel[regLevel].push_back(node_cg);
                                        }

                                    }
                                }
                            }
                        }




                    }

                    // now look for edge neighbors and vertex neighbors of the block, this is only needed when the padding width is >0
                    if(pWidth>0)
                    {
                        const std::vector<unsigned int> blk2Edge_map = blkList[blk].getBlk2DiagMap_vec();
                        const std::vector<unsigned int> blk2Vert_map = blkList[blk].getBlk2VertexMap_vec();
                        const unsigned int blk_ele_1D = blkList[blk].getElemSz1D();

                        for(unsigned int k=0; k< blk_ele_1D; k+=2)
                        {

                            if(blk2Edge_map[2*k]!= LOOK_UP_TABLE_DEFAULT)
                            {
                                if(blk2Edge_map[2*k+0] == blk2Edge_map[2*k+1])
                                {
                                    lookup = blk2Edge_map[2*k+0];
                                    for(unsigned int node = 0; node < nPe; node++)
                                    {
                                        node_cg = e2n_cg[lookup * nPe + node];
                                        if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                            ghostNodesByLevel[regLevel].push_back(node_cg);
                                    }

                                }else
                                {
                                    lookup = blk2Edge_map[2*k+0];
                                    for(unsigned int node = 0; node < nPe; node++)
                                    {
                                        node_cg = e2n_cg[lookup * nPe + node];
                                        if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                            ghostNodesByLevel[regLevel].push_back(node_cg);
                                    }

                                    lookup = blk2Edge_map[2*k+1];
                                    for(unsigned int node = 0; node < nPe; node++)
                                    {
                                        node_cg = e2n_cg[lookup * nPe + node];
                                        if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                            ghostNodesByLevel[regLevel].push_back(node_cg);
                                    }

                                }
                            }
                            
                        }



                        
                        
                        for(unsigned int k=0; k < blk2Vert_map.size(); k++)
                        {
                            lookup = blk2Vert_map[k];
                            if(lookup!=LOOK_UP_TABLE_DEFAULT)
                            {
                                for(unsigned int node = 0; node < nPe; node++)
                                {
                                    node_cg = e2n_cg[lookup * nPe + node];
                                    if(!(node_cg>= nodeLocalBegin && node_cg <nodeLocalEnd))
                                        ghostNodesByLevel[regLevel].push_back(node_cg);
                                }
                            }
                        }



                    }

                }
            
            }

            for(unsigned int l = 0; l < m_uiMaxDepth; l++ )
            {
                std::sort(ghostNodesByLevel[l].begin(), ghostNodesByLevel[l].end());
                ghostNodesByLevel[l].erase( std::unique(ghostNodesByLevel[l].begin(),ghostNodesByLevel[l].end()) , ghostNodesByLevel[l].end() );
            }

            
            const unsigned int * sendSM = &(*(m_uiOctGrid->getSendNodeSM().begin()));
            const unsigned int * recvSM = &(*(m_uiOctGrid->getRecvNodeSM().begin()));

            const unsigned int * sendNodeCounts = &(*(m_uiOctGrid->getNodalSendCounts().begin()));
            const unsigned int * recvNodeCounts = &(*(m_uiOctGrid->getNodalRecvCounts().begin()));

            const unsigned int * sendNodeOffset = &(*(m_uiOctGrid->getNodalSendOffsets().begin()));
            const unsigned int * recvNodeOffset = &(*(m_uiOctGrid->getNodalRecvOffsets().begin()));

            const std::vector<unsigned int> sendNpes = m_uiOctGrid->getSendProcList();
            const std::vector<unsigned int> recvNpes = m_uiOctGrid->getRecvProcList();

            const unsigned int npes = m_uiOctGrid->getMPICommSize();
            
            const unsigned int ghostNodalSz = m_uiOctGrid->getNumPreMeshNodes() + m_uiOctGrid->getNumPostMeshNodes();
            const unsigned int localSz = m_uiOctGrid->getNumLocalMeshNodes();

            m_uiRecvNodeRSM.resize(ghostNodalSz,std::pair<unsigned int,unsigned int>(LOOK_UP_TABLE_DEFAULT,LOOK_UP_TABLE_DEFAULT));
            m_uiSendNodeRSM.resize(localSz,std::pair<unsigned int,unsigned int>(LOOK_UP_TABLE_DEFAULT,LOOK_UP_TABLE_DEFAULT));

            const unsigned int nodePreBegin = m_uiOctGrid->getNodePreGhostBegin();
            const unsigned int nodePreEnd = m_uiOctGrid->getNodePreGhostEnd();

            const unsigned int localBegin = m_uiOctGrid->getNodeLocalBegin();
            const unsigned int localEnd = m_uiOctGrid->getNodeLocalEnd();
            
            const unsigned int nodePostBegin = m_uiOctGrid->getNodePostGhostBegin();
            const unsigned int nodePostEnd = m_uiOctGrid->getNodePostGhostEnd();

            
            for(unsigned int w=0; w< recvNpes.size(); w++)
            {
                const unsigned int p = recvNpes[w];
                for(unsigned int n=recvNodeOffset[p] ; n < (recvNodeOffset[p]+recvNodeCounts[p]); n++)
                {
                    const unsigned int id = recvSM[n];
                    if( id >= nodePreBegin && id < nodePreEnd )
                    {
                        m_uiRecvNodeRSM[(id-nodePreBegin)] = std::pair<unsigned int, unsigned int>(p,n);

                    }else
                    {
                       assert(id>= nodePostBegin && id< nodePostEnd);
                       m_uiRecvNodeRSM[(id-localSz)] = std::pair<unsigned int, unsigned int>(p,n);
                    }
                }
            }

            for(unsigned int w=0; w< sendNpes.size(); w++)
            {
                const unsigned int p = sendNpes[w];
                for(unsigned int n = sendNodeOffset[p] ; n < (sendNodeOffset[p] + sendNodeCounts[p]); n++)
                {
                    const unsigned int id = sendSM[n];
                    m_uiSendNodeRSM[id-nodeLocalBegin] = std::pair<unsigned int, unsigned int>(p,n);
                }
            }

            
            for(unsigned int l=0;l<m_uiMaxDepth;l++)
                ghostNodesByLevel[l].clear();

            delete [] ghostNodesByLevel;




            
            





        }

        
    }

    
}// end of namespace ts
