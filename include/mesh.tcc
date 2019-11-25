//
// Created by milinda on 9/22/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains templated functions in the mesh class.
 * (code moved from mesh.h)
*/
//


namespace ot
{
    template <typename T>
    T* Mesh::createVector() const
    {

        if(!m_uiIsActive) return NULL;

        T* vec=new T[m_uiNumActualNodes];
        return vec;

    }

    template <typename T>
    T* Mesh::createVector(const T initValue) const
    {
        if(!m_uiIsActive) return NULL;

        T* vec=new T[m_uiNumActualNodes];

        for(unsigned int k=0;k<m_uiNumActualNodes;k++)
            vec[k]=initValue;

        return vec;
    }

    template <typename T>
    T* Mesh::createVector(std::function<T(T,T,T)> func) const
    {
        if(!m_uiIsActive) return NULL;

        T* vec=new T[m_uiNumActualNodes];

        unsigned int nodeLookUp_CG;
        unsigned int nodeLookUp_DG;
        unsigned int len;
        double x,y,z;
        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));
        unsigned int ownerID,ii_x,jj_y,kk_z;

        for(unsigned int elem=m_uiElementLocalBegin;elem<m_uiElementLocalEnd;elem++)
        {


            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++ )
                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    {
                        nodeLookUp_CG=m_uiE2NMapping_CG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(nodeLookUp_CG>=m_uiNodeLocalBegin && nodeLookUp_CG<m_uiNodeLocalEnd)
                        {
                            nodeLookUp_DG=m_uiE2NMapping_DG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                            len=1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel());
                            x=pNodes[ownerID].getX()+ ii_x*(len/((double)m_uiElementOrder));
                            y=pNodes[ownerID].getY()+ jj_y*(len/((double)m_uiElementOrder));
                            z=pNodes[ownerID].getZ()+ kk_z*(len/((double)m_uiElementOrder));
                            vec[nodeLookUp_CG]=func(x,y,z);
                        }

                    }

        }

        return vec;
    }

    template <typename T>
    void Mesh::createVector(std::vector<T> &vec) const
    {
        if(!m_uiIsActive) { vec.clear(); return; }
        vec.resize(m_uiNumActualNodes);

    }

    template <typename T>
    void Mesh::createVector(std::vector<T> &vec,const T initValue) const
    {
        if(!m_uiIsActive) { vec.clear(); return; }
        vec.resize(m_uiNumActualNodes,initValue);

    }

    template  <typename T>
    void Mesh::createVector(std::vector<T> & vec, std::function<T(T,T,T)> func) const
    {
        if(!m_uiIsActive) { vec.clear(); return; }
        vec.clear();
        vec.resize(m_uiNumActualNodes,0);
        unsigned int nodeLookUp_CG;
        unsigned int nodeLookUp_DG;
        unsigned int len;
        double x,y,z;
        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));
        unsigned int ownerID,ii_x,jj_y,kk_z;

        for(unsigned int elem=m_uiElementLocalBegin;elem<m_uiElementLocalEnd;elem++)
        {


            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++ )
                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    {
                        nodeLookUp_CG=m_uiE2NMapping_CG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(nodeLookUp_CG>=m_uiNodeLocalBegin && nodeLookUp_CG<m_uiNodeLocalEnd)
                        {
                            nodeLookUp_DG=m_uiE2NMapping_DG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                            len=1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel());
                            x=pNodes[ownerID].getX()+ ii_x*(len/((double)m_uiElementOrder));
                            y=pNodes[ownerID].getY()+ jj_y*(len/((double)m_uiElementOrder));
                            z=pNodes[ownerID].getZ()+ kk_z*(len/((double)m_uiElementOrder));
                            vec[nodeLookUp_CG]=func(x,y,z);
                        }

                    }

        }

    }

    template <typename T>
    void Mesh::createUnZippedVector(std::vector<T> & uvec) const
    {
        if(!m_uiIsActive) { uvec.clear(); return; }
        uvec.resize(m_uiUnZippedVecSz);
    }

    template <typename T>
    void Mesh::createUnZippedVector(std::vector<T> & uvec,const T initValue) const
    {
        if(!m_uiIsActive) { uvec.clear(); return; }
        uvec.resize(m_uiUnZippedVecSz,initValue);
    }


    template <typename T>
    T* Mesh::createUnZippedVector() const
    {
        if(!m_uiIsActive) return NULL;

        T* uvec=new T[m_uiUnZippedVecSz];
        return uvec;
    }

    template <typename T>
    T* Mesh::createUnZippedVector(const T initValue) const
    {
        if(!m_uiIsActive) return NULL;

        T* uvec=new T[m_uiUnZippedVecSz];

        for(unsigned int k=0;k<(m_uiUnZippedVecSz);k++)
            uvec[k]=initValue;

        return uvec;
    }



    template<typename T>
    void Mesh::performGhostExchange(std::vector<T> &vec)
    {

        if((m_uiActiveNpes==1) || (!m_uiIsActive)) return;

        for(unsigned int p=0;p<m_uiActiveNpes;p++) {
            for (unsigned int k = m_uiSendNodeOffset[p]; k < (m_uiSendNodeOffset[p] + m_uiSendNodeCount[p]); k++) {
                m_uiSendBufferNodes[k] = (T)vec[m_uiScatterMapActualNodeSend[k]];
            }
        }

        #ifdef ALLTOALL_SPARSE
            par::Mpi_Alltoallv_sparse(&(*(m_uiSendBufferNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(m_uiRecvBufferNodes.begin())),(int *) (&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),m_uiCommActive);
        #else
            par::Mpi_Alltoallv(&(*(m_uiSendBufferNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(m_uiRecvBufferNodes.begin())),(int *) (&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),m_uiCommActive);
        #endif


        for(unsigned int p=0;p<m_uiActiveNpes;p++)
        {
            for(unsigned int k=m_uiRecvNodeOffset[p];k<(m_uiRecvNodeOffset[p]+m_uiRecvNodeCount[p]);k++)
            {

                //if(fabs(vec[m_uiScatterMapActualNodeRecv[k]]-m_uiRecvBufferNodes[k])>1e-15) std::cout<<"rank: "<<m_uiActiveRank<<" computed: "<<vec[m_uiScatterMapActualNodeRecv[k]]<<" revieved: "<<m_uiRecvBufferNodes[k]<<" recv: from : "<<p<<std::endl;
                vec[m_uiScatterMapActualNodeRecv[k]]=(T)m_uiRecvBufferNodes[k];
            }

        }




    }

    template<typename T>
    void Mesh::performGhostExchange(T*vec)
    {

        if((m_uiActiveNpes==1) || (!m_uiIsActive)) return;

        for(unsigned int p=0;p<m_uiActiveNpes;p++) {
            for (unsigned int k = m_uiSendNodeOffset[p]; k < (m_uiSendNodeOffset[p] + m_uiSendNodeCount[p]); k++) {
                m_uiSendBufferNodes[k] = (T)vec[m_uiScatterMapActualNodeSend[k]];
            }
        }

        #ifdef ALLTOALL_SPARSE
            par::Mpi_Alltoallv_sparse(&(*(m_uiSendBufferNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(m_uiRecvBufferNodes.begin())),(int *) (&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),m_uiCommActive);
        #else
            par::Mpi_Alltoallv(&(*(m_uiSendBufferNodes.begin())),(int *)(&(*(m_uiSendNodeCount.begin()))),(int *)(&(*(m_uiSendNodeOffset.begin()))),&(*(m_uiRecvBufferNodes.begin())),(int *) (&(*(m_uiRecvNodeCount.begin()))),(int *)(&(*(m_uiRecvNodeOffset.begin()))),m_uiCommActive);
        #endif


        for(unsigned int p=0;p<m_uiActiveNpes;p++)
        {
            for(unsigned int k=m_uiRecvNodeOffset[p];k<(m_uiRecvNodeOffset[p]+m_uiRecvNodeCount[p]);k++)
            {

                //if(/*fabs(vec[m_uiScatterMapActualNodeRecv[k]]-m_uiRecvBufferNodes[k])>1e-15*/ isnan(m_uiRecvBufferNodes[k])) std::cout<<"rank: "<<m_uiActiveRank<<" computed: "<<vec[m_uiScatterMapActualNodeRecv[k]]<<" revieved: "<<m_uiRecvBufferNodes[k]<<" recv: from : "<<p<<std::endl;
                vec[m_uiScatterMapActualNodeRecv[k]]=(T)m_uiRecvBufferNodes[k];
            }

        }




    }

    template<typename T>
    void Mesh::ghostExchangeStart(T* vec,T* sendNodeBuffer,T* recvNodeBuffer, MPI_Request * send_reqs, MPI_Request * recv_reqs)
    {

        if((m_uiActiveNpes==1) || (!m_uiIsActive)) return;

        unsigned int proc_id;

        // active recv procs
        for(unsigned int recv_p=0;recv_p<m_uiRecvProcList.size();recv_p++)
        {
            proc_id=m_uiRecvProcList[recv_p];
            recv_reqs[recv_p]=MPI_Request();
            par::Mpi_Irecv((recvNodeBuffer+m_uiRecvNodeOffset[proc_id]),m_uiRecvNodeCount[proc_id],proc_id,0,m_uiCommActive,&recv_reqs[recv_p]);
        }

        for(unsigned int send_p=0;send_p<m_uiSendProcList.size();send_p++)
        {
            proc_id=m_uiSendProcList[send_p];
            for (unsigned int k = m_uiSendNodeOffset[proc_id]; k < (m_uiSendNodeOffset[proc_id] + m_uiSendNodeCount[proc_id]); k++) {
                sendNodeBuffer[k] = (T)vec[m_uiScatterMapActualNodeSend[k]];
            }
        }
        // active send procs
        for(unsigned int send_p=0;send_p<m_uiSendProcList.size();send_p++)
        {
            proc_id=m_uiSendProcList[send_p];
            send_reqs[send_p]=MPI_Request();
            par::Mpi_Isend((sendNodeBuffer+m_uiSendNodeOffset[proc_id]),m_uiSendNodeCount[proc_id],proc_id,0,m_uiCommActive,&send_reqs[send_p]);
        }



    }


    template<typename T>
    void Mesh::ghostExchangeRecvSync(T *vec, T *recvNodeBuffer,MPI_Request *recv_reqs, MPI_Status *recv_sts)
    {
        if((m_uiActiveNpes==1) || (!m_uiIsActive)) return;

        dendro::timer::t_unzip_async_comm.start();
        MPI_Waitall(m_uiRecvProcList.size(),recv_reqs,recv_sts);
        dendro::timer::t_unzip_async_comm.stop();

        unsigned int proc_id=0;
        for(unsigned int recv_p=0;recv_p<m_uiRecvProcList.size();recv_p++)
        {
            proc_id=m_uiRecvProcList[recv_p];
            for(unsigned int k=m_uiRecvNodeOffset[proc_id];k<(m_uiRecvNodeOffset[proc_id]+m_uiRecvNodeCount[proc_id]);k++)
            {
                vec[m_uiScatterMapActualNodeRecv[k]]=(T)recvNodeBuffer[k];
            }
        }

    }


    template <typename T,unsigned int length,unsigned int offsetCentered,unsigned int offsetBackward,unsigned int offsetForward>
    void Mesh::applyStencil(const std::vector<T> &in, std::vector<T> &out,
                            const Stencil<T, length, offsetCentered> &centered,
                            const Stencil<T, length, offsetBackward> &backward,
                            const Stencil<T, length, offsetForward> &forward)
    {


        if(!m_uiIsActive) return;

        double t_uzip;
        double t_uzip_g[3];

        double t_zip;
        double t_zip_g[3];

        double t_stencil;
        double t_stencil_g[3];


        unsigned int blkNpe_1D;
        std::vector<T> unzipVec;
        createUnZippedVector(unzipVec);

        std::vector<T> unzipVec1;
        this->createUnZippedVector(unzipVec1,0.0);

        #ifdef PROFILE_APPLY_STENCIL
                auto t1=std::chrono::high_resolution_clock::now();
        #endif
                this->unzip(&(*(in.begin())),&(*(unzipVec.begin())));
                //std::cout<<"rank: "<<m_uiActiveRank<<" unzip completed "<<std::endl;

        #ifdef PROFILE_APPLY_STENCIL
                auto t2=std::chrono::high_resolution_clock::now();
                t_uzip=std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

                par::Mpi_Reduce(&t_uzip,t_uzip_g,1,MPI_MIN,0,m_uiCommActive);
                par::Mpi_Reduce(&t_uzip,t_uzip_g+1,1,MPI_SUM,0,m_uiCommActive);
                par::Mpi_Reduce(&t_uzip,t_uzip_g+2,1,MPI_MIN,0,m_uiCommActive);
                t_uzip_g[1]=t_uzip_g[1]/(double) m_uiActiveNpes;
        #endif



        unsigned int regLev=0;
        ot::TreeNode blkNode;






        unsigned int centeredOffset=centered.getOffset();
        unsigned int backwardOffset=backward.getOffset();
        unsigned int forwardOffset=forward.getOffset();


        // all the 3 stencil directions should be in the same.
        assert(centered.getStencilDirection()==forward.getStencilDirection());
        assert(centered.getStencilDirection()==backward.getStencilDirection());
        double h=0.0;
        unsigned int lx,ly,lz,offset,paddWidth;
        #ifdef DEBUG_UNZIP_OP
        double d_min=-0.5;
        double d_max=0.5;
        std::function<double(double,double,double)> func =[d_min,d_max](const double x,const double y,const double z){ return (sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
        std::function<double(double,double,double)> dx_func=[d_min,d_max](const double x,const double y,const double z){ return (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(cos(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
        unsigned int  x,y,z,sz,regSz;

        for(unsigned int blk=0;blk<m_uiLocalBlockList.size();blk++) {

            blkNode = m_uiLocalBlockList[blk].getBlockNode();
            regLev = m_uiLocalBlockList[blk].getRegularGridLev();
            lx=m_uiLocalBlockList[blk].getAllocationSzX();
            ly=m_uiLocalBlockList[blk].getAllocationSzY();
            lz=m_uiLocalBlockList[blk].getAllocationSzZ();
            offset=m_uiLocalBlockList[blk].getOffset();
            paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();
            //h=((1u<<(m_uiMaxDepth))*m_uiElementOrder)/((0.5-(-0.5)) * ((1u<<(regLev-blkNode.getLevel()))));
            h = ((blkNode.maxX() - blkNode.minX())) / ((1u << (regLev - blkNode.getLevel())) * m_uiElementOrder);
            h = 1.0 / h;
            blkNpe_1D = m_uiElementOrder * (1u << (regLev - blkNode.getLevel())) + 1 + 2 * paddWidth;
            assert(blkNpe_1D > paddWidth);


            for(unsigned int k=0;k<(blkNpe_1D);k++)
                  for(unsigned int j=0;j<(blkNpe_1D);j++)
                    for(unsigned int i=0;i<(blkNpe_1D);i++)
                    {
                       assert(((blkNode.maxX()-blkNode.minX()))%((1u<<(regLev-blkNode.getLevel()))*m_uiElementOrder)==0);
                       sz=((blkNode.maxX()-blkNode.minX()))/((1u<<(regLev-blkNode.getLevel()))*m_uiElementOrder);
                       regSz=1u<<(m_uiMaxDepth-regLev);

                       if( (i>paddWidth && i<(blkNpe_1D-paddWidth-1)) && (j>paddWidth && j<(blkNpe_1D-paddWidth-1)) && (k>paddWidth && k<(blkNpe_1D-paddWidth-1)))
                       {
                           x=blkNode.getX() + (i-paddWidth)*sz;
                           y=blkNode.getY() + (j-paddWidth)*sz;
                           z=blkNode.getZ() + (k-paddWidth)*sz;
                           if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                               std::cout<<" [internal node mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;

                       }

                        if((blkNode.getX()>=regSz) && (i>=0 && i<=(paddWidth)) && (j>paddWidth && j<(blkNpe_1D-paddWidth)) && (k>paddWidth && k<(blkNpe_1D-paddWidth)))
                        {
                            x=blkNode.getX()-regSz + (i+paddWidth)*sz;
                            y=blkNode.getY() + (j-paddWidth)*sz;
                            z=blkNode.getZ() + (k-paddWidth)*sz;
                            if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                                std::cout<<" [left ghost layer mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;
                        }

                        if((blkNode.getY()>=regSz) && (j>=0 && j<=(paddWidth)) && (i>paddWidth && i<(blkNpe_1D-paddWidth)) && (k>paddWidth && k<(blkNpe_1D-paddWidth)))
                        {
                            x=blkNode.getX() + (i-paddWidth)*sz;
                            y=blkNode.getY()-regSz + (j+paddWidth)*sz;
                            z=blkNode.getZ() + (k-paddWidth)*sz;
                            if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                                std::cout<<" [down ghost layer mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;
                        }


                        if((blkNode.getZ()>=regSz) && (k>=0 && k<=(paddWidth)) && (i>paddWidth && i<(blkNpe_1D-paddWidth)) && (j>paddWidth && j<(blkNpe_1D-paddWidth)))
                        {
                            x=blkNode.getX() + (i-paddWidth)*sz;
                            y=blkNode.getY() + (j-paddWidth)*sz;
                            z=blkNode.getZ()-regSz + (k+paddWidth)*sz;
                            if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                                std::cout<<" [back ghost layer mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;
                        }


                        if((blkNode.maxX()+regSz<=m_uiMeshDomain_max) && (i>=(blkNpe_1D-paddWidth) && i<(blkNpe_1D)) && (j>paddWidth && j<(blkNpe_1D-paddWidth)) && (k>paddWidth && k<(blkNpe_1D-paddWidth)))
                        {
                            x=blkNode.getX()-regSz + (i+paddWidth)*sz;
                            y=blkNode.getY() + (j-paddWidth)*sz;
                            z=blkNode.getZ() + (k-paddWidth)*sz;
                            if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                                std::cout<<" [right ghost layer mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;
                        }


                        if((blkNode.maxY()+regSz<=m_uiMeshDomain_max) && (j>=(blkNpe_1D-paddWidth) && j<(blkNpe_1D)) && (i>paddWidth && i<(blkNpe_1D-paddWidth)) && (k>paddWidth && k<(blkNpe_1D-paddWidth)))
                        {
                            x=blkNode.getX() + (i-paddWidth)*sz;
                            y=blkNode.getY() - regSz + (j+paddWidth)*sz;
                            z=blkNode.getZ() + (k-paddWidth)*sz;
                            if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                                std::cout<<" [up ghost layer mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;
                        }

                        if((blkNode.maxZ()+regSz<=m_uiMeshDomain_max) && (k>=(blkNpe_1D-paddWidth) && k<(blkNpe_1D)) && (i>paddWidth && i<(blkNpe_1D-paddWidth)) && (j>paddWidth && j<(blkNpe_1D-paddWidth)))
                        {
                            x=blkNode.getX() + (i-paddWidth)*sz;
                            y=blkNode.getY() + (j-paddWidth)*sz;
                            z=blkNode.getZ() - regSz + (k+paddWidth)*sz;
                            if(fabs(func(x,y,z)-unzipVec[offset+k*(ly*lx)+j*(lx)+i])>1e-5)
                                std::cout<<" [front ghost layer mismatch] blk: "<<blk<<" blkNode: "<<blkNode<<" sz: "<<sz<<" blkNode_1D: "<<blkNpe_1D<<" (x,y,z): ( "<<x<<", "<<y<<", "<<z<<")  (i,j,k)= ("<<i<<","<<j<<", "<<k<<" )"<<" ) func: "<<func(x,y,z)<<" : read value "<<unzipVec[offset+k*(ly*lx)+j*(lx)+i]<<std::endl;
                        }




                    }



        }

        #endif

        #ifdef PROFILE_APPLY_STENCIL
                t1=std::chrono::high_resolution_clock::now();
        #endif



        if(centered.getStencilDirection()==StencilDirection::STENCIL_DIR_X)
        {

            for(unsigned int blk=0;blk<m_uiLocalBlockList.size();blk++)
            {
                blkNode=m_uiLocalBlockList[blk].getBlockNode();
                regLev=m_uiLocalBlockList[blk].getRegularGridLev();

                lx=m_uiLocalBlockList[blk].getAllocationSzX();
                ly=m_uiLocalBlockList[blk].getAllocationSzY();
                lz=m_uiLocalBlockList[blk].getAllocationSzZ();
                offset=m_uiLocalBlockList[blk].getOffset();
                paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();

                //h=((1u<<(m_uiMaxDepth))*m_uiElementOrder)/((0.5-(-0.5)) * ((1u<<(regLev-blkNode.getLevel()))));
                h=((blkNode.maxX()-blkNode.minX()))/((double)(1u<<(regLev-blkNode.getLevel()))*m_uiElementOrder);
                h=1.0/h;
                blkNpe_1D=m_uiElementOrder*(1u<<(regLev-blkNode.getLevel()))+1+2*paddWidth;
                assert(blkNpe_1D>paddWidth);


                if(blkNode.minX()==m_uiMeshDomain_min)
                {

                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying forward difference difference: "<<std::endl;
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=paddWidth;i<2*paddWidth;i++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<forward.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=forward[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+i+index-forwardOffset]*h;

                            }


                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=2*paddWidth;i<(blkNpe_1D-2*paddWidth);i++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+i+index-centeredOffset]*h;

                            }

                    if(blkNode.maxX()==m_uiMeshDomain_max)
                    {
                        for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                            for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                                for(unsigned int i=(blkNpe_1D-2*paddWidth);i<(blkNpe_1D-paddWidth);i++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                    for(unsigned int index=0;index<backward.getStencilLength();index++)
                                        unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=backward[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+i+index-backwardOffset]*h;
                                }

                    }else
                    {
                        assert(blkNode.maxX()<m_uiMeshDomain_max);
                        for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                            for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                                for(unsigned int i=(blkNpe_1D-2*paddWidth);i<(blkNpe_1D-paddWidth);i++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                    for(unsigned int index=0;index<centered.getStencilLength();index++)
                                        unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+i+index-centeredOffset]*h;

                                }
                    }



                }else if(blkNode.maxX()==m_uiMeshDomain_max)
                {

                    assert(blkNode.minX()>m_uiMeshDomain_min);
                    assert((blkNpe_1D-2*paddWidth));
                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying backward difference difference: "<<std::endl;
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=paddWidth;i<(blkNpe_1D-2*paddWidth);i++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+i+index-centeredOffset]*h;

                            }


                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=(blkNpe_1D-2*paddWidth);i<(blkNpe_1D-paddWidth);i++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<backward.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=backward[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+i+index-backwardOffset]*h;
                            }


                }else
                {
                    assert(blkNode.minX()>m_uiMeshDomain_min && blkNode.maxX()<m_uiMeshDomain_max);
                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying centered difference difference: "<<std::endl;
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+j*(lx)+(i+index-centeredOffset)]*h;

                                }

                            }

                }


            }

        }


        if(centered.getStencilDirection()==StencilDirection::STENCIL_DIR_Y)
        {

            for(unsigned int blk=0;blk<m_uiLocalBlockList.size();blk++)
            {
                blkNode=m_uiLocalBlockList[blk].getBlockNode();
                regLev=m_uiLocalBlockList[blk].getRegularGridLev();
                //h=((1u<<(m_uiMaxDepth))*m_uiElementOrder)/((0.5-(-0.5)) * ((1u<<(regLev-blkNode.getLevel()))));

                lx=m_uiLocalBlockList[blk].getAllocationSzX();
                ly=m_uiLocalBlockList[blk].getAllocationSzY();
                lz=m_uiLocalBlockList[blk].getAllocationSzZ();
                offset=m_uiLocalBlockList[blk].getOffset();
                paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();


                h=((blkNode.maxY()-blkNode.minY()))/((double)(1u<<(regLev-blkNode.getLevel()))*m_uiElementOrder);
                h=1.0/h;
                blkNpe_1D=m_uiElementOrder*(1u<<(regLev-blkNode.getLevel()))+1+2*paddWidth;
                assert(blkNpe_1D>paddWidth);


                if(blkNode.minY()==m_uiMeshDomain_min)
                {

                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying forward difference difference: "<<std::endl;
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int j=paddWidth;j<2*paddWidth;j++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<forward.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=forward[index]*unzipVec[offset+k*(ly*lx)+(j+index-forwardOffset)*(lx)+i]*h;

                            }


                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int j=2*paddWidth;j<(blkNpe_1D-2*paddWidth);j++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+(j+index-centeredOffset)*(lx)+i]*h;

                            }

                    if(blkNode.maxY()==m_uiMeshDomain_max)
                    {
                        for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                            for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                                for(unsigned int j=(blkNpe_1D-2*paddWidth);j<(blkNpe_1D-paddWidth);j++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                    for(unsigned int index=0;index<backward.getStencilLength();index++)
                                        unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=backward[index]*unzipVec[offset+k*(ly*lx)+(j+index-backwardOffset)*(lx)+i]*h;
                                }

                    }else
                    {
                        assert(blkNode.maxY()<m_uiMeshDomain_max);
                        for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                            for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                                for(unsigned int j=(blkNpe_1D-2*paddWidth);j<(blkNpe_1D-paddWidth);j++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                    for(unsigned int index=0;index<centered.getStencilLength();index++)
                                        unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+(j+index-centeredOffset)*(lx)+i]*h;

                                }
                    }



                }else if(blkNode.maxY()==m_uiMeshDomain_max)
                {

                    assert(blkNode.minY()>m_uiMeshDomain_min);
                    assert((blkNpe_1D-2*paddWidth));
                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying backward difference difference: "<<std::endl;
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int j=paddWidth;j<(blkNpe_1D-2*paddWidth);j++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+(j+index-centeredOffset)*(lx)+i]*h;

                            }


                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int j=(blkNpe_1D-2*paddWidth);j<(blkNpe_1D-paddWidth);j++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<backward.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=backward[index]*unzipVec[offset+k*(ly*lx)+(j+index-backwardOffset)*(lx)+i]*h;
                            }


                }else
                {
                    assert(blkNode.minY()>m_uiMeshDomain_min && blkNode.maxY()<m_uiMeshDomain_max);
                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying centered difference difference: "<<std::endl;
                    for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+k*(ly*lx)+(j+index-centeredOffset)*(lx)+(i)]*h;

                                }

                            }

                }


            }

        }




        if(centered.getStencilDirection()==StencilDirection::STENCIL_DIR_Z)
        {


            for(unsigned int blk=0;blk<m_uiLocalBlockList.size();blk++)
            {
                blkNode=m_uiLocalBlockList[blk].getBlockNode();
                regLev=m_uiLocalBlockList[blk].getRegularGridLev();
                //h=((1u<<(m_uiMaxDepth))*m_uiElementOrder)/((0.5-(-0.5)) * ((1u<<(regLev-blkNode.getLevel()))));

                lx=m_uiLocalBlockList[blk].getAllocationSzX();
                ly=m_uiLocalBlockList[blk].getAllocationSzY();
                lz=m_uiLocalBlockList[blk].getAllocationSzZ();
                offset=m_uiLocalBlockList[blk].getOffset();
                paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();

                h=((blkNode.maxZ()-blkNode.minZ()))/((double)(1u<<(regLev-blkNode.getLevel()))*m_uiElementOrder);
                h=1.0/h;
                blkNpe_1D=m_uiElementOrder*(1u<<(regLev-blkNode.getLevel()))+1+2*paddWidth;
                assert(blkNpe_1D>paddWidth);


                if(blkNode.minZ()==m_uiMeshDomain_min)
                {

                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying forward difference difference: "<<std::endl;
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int k=paddWidth;k<2*paddWidth;k++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<forward.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=forward[index]*unzipVec[offset+(k+index-forwardOffset)*(ly*lx)+(j)*(lx)+i]*h;

                            }


                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int k=2*paddWidth;k<(blkNpe_1D-2*paddWidth);k++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+(k+index-centeredOffset)*(ly*lx)+(j)*(lx)+i]*h;

                            }

                    if(blkNode.maxZ()==m_uiMeshDomain_max)
                    {
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                                for(unsigned int k=(blkNpe_1D-2*paddWidth);k<(blkNpe_1D-paddWidth);k++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                    for(unsigned int index=0;index<backward.getStencilLength();index++)
                                        unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=backward[index]*unzipVec[offset+(k+index-backwardOffset)*(ly*lx)+(j)*(lx)+i]*h;
                                }

                    }else
                    {
                        assert(blkNode.maxZ()<m_uiMeshDomain_max);
                        for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                            for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                                for(unsigned int k=(blkNpe_1D-2*paddWidth);k<(blkNpe_1D-paddWidth);k++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                    for(unsigned int index=0;index<centered.getStencilLength();index++)
                                        unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+(k+index-centeredOffset)*(ly*lx)+(j)*(lx)+i]*h;

                                }
                    }



                }else if(blkNode.maxZ()==m_uiMeshDomain_max)
                {

                    assert(blkNode.minZ()>m_uiMeshDomain_min);
                    assert((blkNpe_1D-2*paddWidth));
                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying backward difference difference: "<<std::endl;
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int k=paddWidth;k<(blkNpe_1D-2*paddWidth);k++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+(k+index-centeredOffset)*(ly*lx)+(j)*(lx)+i]*h;

                            }


                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int k=(blkNpe_1D-2*paddWidth);k<(blkNpe_1D-paddWidth);k++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<backward.getStencilLength();index++)
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=backward[index]*unzipVec[offset+(k+index-backwardOffset)*(ly*lx)+(j)*(lx)+i]*h;
                            }


                }else
                {
                    assert(blkNode.minZ()>m_uiMeshDomain_min && blkNode.maxZ()<m_uiMeshDomain_max);
                    //std::cout<<"rank: "<<m_uiActiveRank<<" applying centered difference difference: "<<std::endl;
                    for(unsigned int j=paddWidth;j<(blkNpe_1D-paddWidth);j++)
                        for(unsigned int i=paddWidth;i<(blkNpe_1D-paddWidth);i++)
                            for(unsigned int k=paddWidth;k<(blkNpe_1D-paddWidth);k++)
                            {
                                unzipVec1[offset+k*(ly*lx)+j*(lx)+i]=0;
                                for(unsigned int index=0;index<centered.getStencilLength();index++)
                                {
                                    unzipVec1[offset+k*(ly*lx)+j*(lx)+i]+=centered[index]*unzipVec[offset+(k+index-centeredOffset)*(ly*lx)+(j)*(lx)+(i)]*h;

                                }

                            }

                }


            }


        }

        #ifdef PROFILE_APPLY_STENCIL
                t2=std::chrono::high_resolution_clock::now();
                t_stencil=std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

                par::Mpi_Reduce(&t_stencil,t_stencil_g,1,MPI_MIN,0,m_uiCommActive);
                par::Mpi_Reduce(&t_stencil,t_stencil_g+1,1,MPI_SUM,0,m_uiCommActive);
                par::Mpi_Reduce(&t_stencil,t_stencil_g+2,1,MPI_MIN,0,m_uiCommActive);
                t_stencil_g[1]=t_stencil_g[1]/(double) m_uiActiveNpes;

                t1=std::chrono::high_resolution_clock::now();
        #endif
                this->createVector(out);
                this->zip(&(*(unzipVec1.begin())),&(*(out.begin())));

        #ifdef PROFILE_APPLY_STENCIL
                t2=std::chrono::high_resolution_clock::now();
                t_zip=std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
                par::Mpi_Reduce(&t_zip,t_zip_g,1,MPI_MIN,0,m_uiCommActive);
                par::Mpi_Reduce(&t_zip,t_zip_g+1,1,MPI_SUM,0,m_uiCommActive);
                par::Mpi_Reduce(&t_zip,t_zip_g+2,1,MPI_MIN,0,m_uiCommActive);
                t_zip_g[1]=t_zip_g[1]/(double) m_uiActiveNpes;

                if(!m_uiActiveRank)
                {
                    std::cout<<"unzip_max \t stencil_max \t zip_max "<<std::endl;
                    std::cout<<t_uzip_g[1]<<" \t "<<t_stencil_g[1]<<" \t "<<t_zip_g[1]<<std::endl;
                }

        #endif
        unzipVec1.clear();
        unzipVec.clear();

    }





    template<typename T>
    inline void Mesh::interpDownWind(const double *downWind,const unsigned int element,const unsigned int lookup, T* vecLookUp, const unsigned int cnum, const T* parentInterpIn, T* parentInterpOut,const unsigned int padDir,const unsigned int padWidth,const T* zippedVec,T* out)
    {

        if(!m_uiIsActive) return;

        // no need for interpolation if the padding width is less than eq. 2
        if(padWidth<=2) return;
        if(m_uiElementOrder!=4) return; // written only for order 4 elements.

        this->parent2ChildInterpolation(parentInterpIn,parentInterpOut,cnum,3);
        this->getElementNodalValues(zippedVec,vecLookUp,lookup);

        T interpVal;

        #ifdef DEBUG_UPWIND_INTERP
                #pragma message("DEBIG_DOWNWIND_INTERP: ON")
                T adv_val1,adv_val2;
        #endif
        const unsigned int stencilWidth=5;
        if(padDir==OCT_DIR_LEFT)
        {
            for(unsigned int k=0;k<(m_uiElementOrder+1);k+=2)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j+=2)
                {
                    interpVal=2*downWind[0]*vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+(1)];
                    for(unsigned int index=0;index<(stencilWidth-2);index++)
                        interpVal+=((downWind[index+1])*(2*vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+(index+2)]-vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+(2*index)]));

                    interpVal+=((downWind[4])*(2*parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+1]-parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+2]));
                    interpVal/=downWind[0];

                    out[((((cnum & 4u)>>2u)*m_uiElementOrder+k)>>1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((((cnum & 2u)>>1u)*m_uiElementOrder+j)>>1)*(m_uiElementOrder+1)+(((((cnum & (1u)))*m_uiElementOrder+2)>>1))]=interpVal;

            #ifdef DEBUG_UPWIND_INTERP
                    adv_val1=2*downWind[4]*parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+1];
                    for(unsigned int index=0;index<(stencilWidth-1);index++)
                        adv_val1+=2*vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+index+1]*downWind[index];

                    adv_val2=downWind[4]*parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+2];
                    adv_val2+=downWind[0]*interpVal;
                    for(unsigned int index=1;index<(stencilWidth-1);index++)
                        adv_val2+=vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+(2*index-2)]*downWind[index];

                    if(fabs(adv_val1-adv_val2)>1e-3)
                        std::cout<<"[left] m_uiActiveRank: "<<m_uiActiveRank<<" adv_val1: "<<adv_val1<<" adv_val_2: "<<adv_val2<<std::endl;

            #endif

                }



        }else if (padDir==OCT_DIR_DOWN)
        {

            for(unsigned int k=0;k<(m_uiElementOrder+1);k+=2)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i+=2)
                {
                    interpVal=2*downWind[0]*vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(1)*(m_uiElementOrder+1)+(i)];
                    for(unsigned int index=0;index<(stencilWidth-2);index++)
                        interpVal+=((downWind[index+1])*(2*vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(index+2)*(m_uiElementOrder+1)+(i)]-vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(2*index)*(m_uiElementOrder+1)+i]));

                    interpVal+=((downWind[4])*(2*parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(1)*(m_uiElementOrder+1)+i]-parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(2)*(m_uiElementOrder+1)+i]));
                    interpVal/=downWind[0];

                    out[((((cnum & 4u)>>2u)*m_uiElementOrder+k)>>1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((((cnum & 2u)>>1u)*m_uiElementOrder+2)>>1)*(m_uiElementOrder+1)+(((((cnum & (1u)))*m_uiElementOrder+i)>>1))]=interpVal;

                #ifdef DEBUG_UPWIND_INTERP
                    adv_val1=2*downWind[4]*parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(1)*(m_uiElementOrder+1)+i];
                    for(unsigned int index=0;index<(stencilWidth-1);index++)
                        adv_val1+=2*vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(index+1)*(m_uiElementOrder+1)+i]*downWind[index];

                    adv_val2=downWind[4]*parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+2*(m_uiElementOrder+1)+i];
                    adv_val2+=downWind[0]*interpVal;
                    for(unsigned int index=1;index<(stencilWidth-1);index++)
                        adv_val2+=vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(2*index-2)*(m_uiElementOrder+1)+i]*downWind[index];

                    if(fabs(adv_val1-adv_val2)>1e-3)
                        std::cout<<"[down] m_uiActiveRank: "<<m_uiActiveRank<<" adv_val1: "<<adv_val1<<" adv_val_2: "<<adv_val2<<std::endl;

                #endif
                }



        }else if (padDir==OCT_DIR_BACK)
        {

            for(unsigned int j=0;j<(m_uiElementOrder+1);j+=2)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i+=2)
                {
                    interpVal=2*downWind[0]*vecLookUp[(1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+(i)];
                    for(unsigned int index=0;index<(stencilWidth-2);index++)
                        interpVal+=((downWind[index+1])*(2*vecLookUp[(index+2)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+(i)]-vecLookUp[(2*index)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]));

                    interpVal+=((downWind[4])*(2*parentInterpOut[(1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]-parentInterpOut[(2)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]));
                    interpVal/=downWind[0];

                    out[((((cnum & 4u)>>2u)*m_uiElementOrder+2)>>1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((((cnum & 2u)>>1u)*m_uiElementOrder+j)>>1)*(m_uiElementOrder+1)+(((((cnum & (1u)))*m_uiElementOrder+i)>>1))]=interpVal;

                #ifdef DEBUG_UPWIND_INTERP
                    adv_val1=2*downWind[4]*parentInterpOut[(1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    for(unsigned int index=0;index<(stencilWidth-1);index++)
                        adv_val1+=2*vecLookUp[(index+1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]*downWind[index];

                    adv_val2=downWind[4]*parentInterpOut[(2)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    adv_val2+=downWind[0]*interpVal;
                    for(unsigned int index=1;index<(stencilWidth-1);index++)
                        adv_val2+=vecLookUp[(2*index-2)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]*downWind[index];

                    if(fabs(adv_val1-adv_val2)>1e-3)
                        std::cout<<"[back] m_uiActiveRank: "<<m_uiActiveRank<<" adv_val1: "<<adv_val1<<" adv_val_2: "<<adv_val2<<std::endl;

                #endif
                }


        }else
        {
            std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" invalid padding direction specified "<<std::endl;
        }


    }



    template<typename T>
    inline void Mesh::interpUpWind(const double *upWind,const unsigned int element,const unsigned int lookup, T* vecLookUp, const unsigned int cnum, const T* parentInterpIn, T* parentInterpOut,const unsigned int padDir,const unsigned int padWidth,const T* zippedVec,T* out)
    {

        if(!m_uiIsActive) return;

        // no need for interpolation if the padding width is less than eq. 2
        if(padWidth<=2) return;
        if(m_uiElementOrder!=4) return; // written only for order 4 elements.

        this->parent2ChildInterpolation(parentInterpIn,parentInterpOut,cnum,3);
        this->getElementNodalValues(zippedVec,vecLookUp,lookup);

        const unsigned int stencilWidth=5;
        T interpVal;

        #ifdef DEBUG_UPWIND_INTERP
        #pragma message("DEBIG_UPWIND_INTERP: ON")
                T adv_val1;
                T adv_val2;
        #endif


        if(padDir==OCT_DIR_RIGHT)
        {
            for(unsigned int k=0;k<(m_uiElementOrder+1);k+=2)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j+=2)
                {

                    interpVal=2*upWind[4]*vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+(3)];
                    interpVal+=((upWind[0])*(2*parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+3]-parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+2]));
                    for(unsigned int index=0;index<(stencilWidth-2);index++)
                        interpVal+=((upWind[index+1])*(2*vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+index]-vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+(2*index)]));

                    interpVal/=upWind[4];
                    out[((((cnum & 4u)>>2u)*m_uiElementOrder+k)>>1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((((cnum & 2u)>>1u)*m_uiElementOrder+j)>>1)*(m_uiElementOrder+1)+(((((cnum & (1u)))*m_uiElementOrder+2)>>1))]=interpVal;

                #ifdef DEBUG_UPWIND_INTERP
                    adv_val1=2*upWind[0]*parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+3];
                    for(unsigned int index=0;index<(stencilWidth-1);index++)
                        adv_val1+=2*vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+index]*upWind[index+1];

                    adv_val2=upWind[0]*parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+2];
                    for(unsigned int index=1;index<(stencilWidth-1);index++)
                        adv_val2+=vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+(2*index-2)]*upWind[index];

                    adv_val2+=upWind[4]*interpVal;

                    if(fabs(adv_val1-adv_val2)>1e-3)
                        std::cout<<"[right] m_uiActiveRank: "<<m_uiActiveRank<<" adv_val1: "<<adv_val1<<" adv_val_2: "<<adv_val2<<std::endl;

                #endif

                }






        }else if (padDir==OCT_DIR_UP)
        {

            for(unsigned int k=0;k<(m_uiElementOrder+1);k+=2)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i+=2)
                {
                    interpVal=2*upWind[4]*vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(3)*(m_uiElementOrder+1)+i];
                    interpVal+=((upWind[0])*(2*parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(3)*(m_uiElementOrder+1)+i]-parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(2)*(m_uiElementOrder+1)+i]));
                    for(unsigned int index=0;index<(stencilWidth-2);index++)
                        interpVal+=((upWind[index+1])*(2*vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(index)*(m_uiElementOrder+1)+i]-vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(2*index)*(m_uiElementOrder+1)+i]));

                    interpVal/=upWind[4];
                    out[((((cnum & 4u)>>2u)*m_uiElementOrder+k)>>1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((((cnum & 2u)>>1u)*m_uiElementOrder+2)>>1)*(m_uiElementOrder+1)+(((((cnum & (1u)))*m_uiElementOrder+i)>>1))]=interpVal;

                #ifdef DEBUG_UPWIND_INTERP
                    adv_val1=2*upWind[0]*parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(3)*(m_uiElementOrder+1)+i];
                    for(unsigned int index=0;index<(stencilWidth-1);index++)
                        adv_val1+=2*vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(index)*(m_uiElementOrder+1)+i]*upWind[index+1];

                    adv_val2=upWind[0]*parentInterpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(2)*(m_uiElementOrder+1)+i];
                    for(unsigned int index=1;index<(stencilWidth-1);index++)
                        adv_val2+=vecLookUp[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+(2*index-2)*(m_uiElementOrder+1)+i]*upWind[index];

                    adv_val2+=upWind[4]*interpVal;

                    if(fabs(adv_val1-adv_val2)>1e-3)
                        std::cout<<"[up] m_uiActiveRank: "<<m_uiActiveRank<<" adv_val1: "<<adv_val1<<" adv_val_2: "<<adv_val2<<std::endl;

                #endif

                }


        }else if (padDir==OCT_DIR_FRONT)
        {

            for(unsigned int j=0;j<(m_uiElementOrder+1);j+=2)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i+=2)
                {
                    interpVal=2*upWind[4]*vecLookUp[(3)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    interpVal+=((upWind[0])*(2*parentInterpOut[(3)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]-parentInterpOut[(2)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]));
                    for(unsigned int index=0;index<(stencilWidth-2);index++)
                        interpVal+=((upWind[index+1])*(2*vecLookUp[(index)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]-vecLookUp[(2*index)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]));

                    interpVal/=upWind[4];
                    out[((((cnum & 4u)>>2u)*m_uiElementOrder+2)>>1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((((cnum & 2u)>>1u)*m_uiElementOrder+j)>>1)*(m_uiElementOrder+1)+(((((cnum & (1u)))*m_uiElementOrder+i)>>1))]=interpVal;

                #ifdef DEBUG_UPWIND_INTERP
                    adv_val1=2*upWind[0]*parentInterpOut[(3)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    for(unsigned int index=0;index<(stencilWidth-1);index++)
                        adv_val1+=2*vecLookUp[(index)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]*upWind[index+1];

                    adv_val2=upWind[0]*parentInterpOut[(2)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    for(unsigned int index=1;index<(stencilWidth-1);index++)
                        adv_val2+=vecLookUp[(2*index-2)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]*upWind[index];

                    adv_val2+=upWind[4]*interpVal;

                    if(fabs(adv_val1-adv_val2)>1e-3)
                        std::cout<<"[front] m_uiActiveRank: "<<m_uiActiveRank<<" adv_val1: "<<adv_val1<<" adv_val_2: "<<adv_val2<<std::endl;

                #endif
                }

        }else
        {
            std::cout<<"m_uiActiveRank: "<<m_uiActiveRank<<" invalid padding direction specified "<<std::endl;
        }


    }


    template <typename pKey, typename pNode>
    void Mesh::searchKeys(std::vector<pKey>& pKeys,std::vector<pNode>& pNodes)
    {


        assert(seq::test::isSorted(pNodes));

        std::vector<Key> pKeys_cpy;
        pKeys_cpy.resize(pKeys.size());

        for(unsigned int k=0;k<pKeys.size();k++)
        {
            pKeys_cpy[k]=pKeys[k];
            pKeys_cpy[k].addOwner(k);
            pKeys_cpy[k].setSearchResult(LOOK_UP_TABLE_DEFAULT);
        }

        SFC::seqSearch::SFC_treeSearch(&(*(pKeys_cpy.begin())),&(*(pNodes.begin())),0,pKeys_cpy.size(),0,pNodes.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);

        for(unsigned int k=0;k<pKeys_cpy.size();k++)
        {
            if((pKeys_cpy[k].getFlag() & OCT_FOUND)) {
                pKeys[(*(pKeys_cpy[k].getOwnerList()))[0]].setSearchResult(pKeys_cpy[k].getSearchResult());
                pKeys[k].setFlag((pKeys[k].getFlag() | OCT_FOUND));
            }

        }

        pKeys_cpy.clear();




    }


   
    template <typename T>
    bool Mesh::isReMeshUnzip(const T **unzippedVec,const unsigned int * varIds,const unsigned int numVars,std::function<double(double,double,double)>wavelet_tol,double amr_coarse_fac, double coarsen_hx)
    {

        bool isOctChange=false;
        if(m_uiIsActive)
        {
            // remove all the previously set falgs if there is any.  THIS will change the all flags to no CHANGE
            for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
                m_uiAllElements[ele].setFlag(((OCT_NO_CHANGE<<NUM_LEVEL_BITS)|m_uiAllElements[ele].getLevel()));

            ot::TreeNode blkNode;
            unsigned int sz[3];
            double dh[3];
            unsigned int bflag,offset;
            unsigned int regLev;
            unsigned int eIndex[3];
            double *  waveletR = NULL;
            double *  waveletC = NULL;
            unsigned int num_wr =0 ,num_wc =0;

            double * wsIn = new double[m_uiNpE];
            double * wsOut = new double[m_uiNpE];
            double ** ws = new double*[2];
            ws[0] = wsIn;
            ws[1] = wsOut;

            // upper bound for the refine and coarsen wavelets.     
            waveletR = new double[64];
            num_wr = 64;
            
            waveletC = new double[64];
            num_wc = 64 ;

            const unsigned int paddWidth=3;
            unsigned int eleIndexMin=0,eleIndexMax=0;

            double l_inf;
            double x,y,z,tol;


            // first pass to identify the refined elements.
            for(unsigned blk=0;blk<m_uiLocalBlockList.size();blk++)
            {

                blkNode=m_uiLocalBlockList[blk].getBlockNode();

                sz[0]=m_uiLocalBlockList[blk].getAllocationSzX();
                sz[1]=m_uiLocalBlockList[blk].getAllocationSzY();
                sz[2]=m_uiLocalBlockList[blk].getAllocationSzZ();

                bflag=m_uiLocalBlockList[blk].getBlkNodeFlag();
                offset=m_uiLocalBlockList[blk].getOffset();

                regLev=m_uiLocalBlockList[blk].getRegularGridLev();
                eleIndexMax=(1u<<(regLev-blkNode.getLevel()))-1;

                //if(bflag!=0) continue;

                for(unsigned int ele=m_uiLocalBlockList[blk].getLocalElementBegin();ele<m_uiLocalBlockList[blk].getLocalElementEnd();ele++)
                {

                    if((m_uiAllElements[ele].getLevel()+MAXDEAPTH_LEVEL_DIFF+1)>=m_uiMaxDepth) continue;

                    x=m_uiAllElements[ele].getX();
                    y=m_uiAllElements[ele].getY();
                    z=m_uiAllElements[ele].getZ();
                    tol=wavelet_tol(x,y,z);

                    eIndex[0]=(m_uiAllElements[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                    eIndex[1]=(m_uiAllElements[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                    eIndex[2]=(m_uiAllElements[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                    if((bflag &(1u<<OCT_DIR_LEFT)) && eIndex[0]==eleIndexMin)   continue;
                    if((bflag &(1u<<OCT_DIR_DOWN)) && eIndex[1]==eleIndexMin)   continue;
                    if((bflag &(1u<<OCT_DIR_BACK)) && eIndex[2]==eleIndexMin)   continue;

                    if((bflag &(1u<<OCT_DIR_RIGHT)) && eIndex[0]==eleIndexMax)  continue;
                    if((bflag &(1u<<OCT_DIR_UP)) && eIndex[1]==eleIndexMax)     continue;
                    if((bflag &(1u<<OCT_DIR_FRONT)) && eIndex[2]==eleIndexMax)  continue;

                    for(unsigned int var=0;var<numVars;var++)
                    {

                         refine_wavelets(&unzippedVec[varIds[var]][offset],m_uiElementOrder,eIndex,paddWidth,sz,waveletR,num_wr,(double**)ws);
                        //  for(unsigned int k=0; k<4; k+=3)
                        //      for(unsigned int j=0; j<4; j+=3)
                        //       for(unsigned int i=0; i<4; i+=3)                                
                        //         waveletR[k*16 + j*4 + i] =0;

                         //l_inf=normLInfty(waveletR,num_wr);
                         l_inf = normL2(waveletR,num_wr)/num_wr;

                            // for(unsigned int k=1; k<3; k+=1)
                            //   for(unsigned int j=1; j<3; j+=1)
                            //    for(unsigned int i=1; i<3; i+=1)
                            //     std::cout<<"ref1: (i,j,k) : " << (i-1)<<" , "<<(j-1)<<" , "<<(k-1)<<": "<<waveletR[k*16 + j*4 + i]<<std::endl;
                            
                        
                       
                        // computeRefineWavelets(&unzippedVec[varIds[var]][offset],0,m_uiElementOrder,eIndex,paddWidth,sz,waveletR);
                        // l_inf=normLInfty(waveletR,NUM_REFINE_WAVELET_COEF);

                        //     for(unsigned int k=1; k<3; k+=1)
                        //       for(unsigned int j=1; j<3; j+=1)
                        //        for(unsigned int i=1; i<3; i+=1)
                        //         std::cout<<"ref2: (i,j,k) : " << (i-1)<<" , "<<(j-1)<<" , "<<(k-1)<<": "<<waveletR[(k-1)*4 + (j-1)*2 +i-1]<<std::endl;

                        if(l_inf>tol)
                        {
                            // for(unsigned int k=0;k<num_wr;k++)
                            //    std::cout<<"elem: "<<m_uiAllElements[ele]<<" wr["<<k<<"]: "<<waveletR[k]<<std::endl;
                            assert((m_uiAllElements[ele].getLevel()+MAXDEAPTH_LEVEL_DIFF+1)<m_uiMaxDepth);
                            //std::cout<<"rank: "<<m_uiActiveRank<<" element R: "<<m_uiAllElements[ele]<<" w_tol: "<<l_inf<<std::endl;
                            m_uiAllElements[ele].setFlag(((OCT_SPLIT<<NUM_LEVEL_BITS)|m_uiAllElements[ele].getLevel()));
                            assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT);
                            break; // no point checking for the other variables since this octree needs to be refined.
                        }

                    }



                }


            }

            //second pass to identify the coarsening elements.
            for(unsigned blk=0;blk<m_uiLocalBlockList.size();blk++)
            {

                blkNode=m_uiLocalBlockList[blk].getBlockNode();

                sz[0]=m_uiLocalBlockList[blk].getAllocationSzX();
                sz[1]=m_uiLocalBlockList[blk].getAllocationSzY();
                sz[2]=m_uiLocalBlockList[blk].getAllocationSzZ();

                bflag=m_uiLocalBlockList[blk].getBlkNodeFlag();
                offset=m_uiLocalBlockList[blk].getOffset();

                regLev=m_uiLocalBlockList[blk].getRegularGridLev();
                eleIndexMax=(1u<<(regLev-blkNode.getLevel()))-1;

                dh[0]=coarsen_hx*(m_uiLocalBlockList[blk].computeGridDx());
                dh[1]=coarsen_hx*(m_uiLocalBlockList[blk].computeGridDy());
                dh[2]=coarsen_hx*(m_uiLocalBlockList[blk].computeGridDz());


                if((eleIndexMax==0) || (bflag!=0)) continue; // this implies the blocks with only 1 child and boundary blocks.

                bool isEligibleCoarsen=true;
                bool isCoarsen=true;
                ot::TreeNode tmpOct;

                for(unsigned int ele=m_uiLocalBlockList[blk].getLocalElementBegin();ele<m_uiLocalBlockList[blk].getLocalElementEnd();ele+=NUM_CHILDREN)
                {

                    assert(m_uiAllElements[ele].getParent()==m_uiAllElements[ele+NUM_CHILDREN-1].getParent());

                    isEligibleCoarsen=true;
                    for(unsigned int child=0;child<NUM_CHILDREN;child++)
                    {
                        if((m_uiAllElements[ele+child].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                        {
                            isEligibleCoarsen=false;
                            break;
                        }

                    }

                    if((isEligibleCoarsen) && (m_uiAllElements[ele].getLevel()>1))
                    {
                        tmpOct=m_uiAllElements[ele].getParent();
                        x=tmpOct.getX() + (1u<<(m_uiMaxDepth-tmpOct.getLevel()-1));
                        y=tmpOct.getY() + (1u<<(m_uiMaxDepth-tmpOct.getLevel()-1));
                        z=tmpOct.getZ() + (1u<<(m_uiMaxDepth-tmpOct.getLevel()-1));
                        tol=wavelet_tol(x,y,z);
                        tmpOct=ot::TreeNode(tmpOct.getX(),tmpOct.getY(),tmpOct.getZ(),tmpOct.getLevel()+1,m_uiDim,m_uiMaxDepth);

                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                        {
                            if(tmpOct==m_uiAllElements[ele+child])
                            {
                                eIndex[0]=(m_uiAllElements[ele+child].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                                eIndex[1]=(m_uiAllElements[ele+child].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                                eIndex[2]=(m_uiAllElements[ele+child].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                                break;
                            }

                        }

                        isCoarsen=true;

                        for(unsigned int var=0;var<numVars;var++)
                        {
                            coarsen_wavelets(&unzippedVec[varIds[var]][offset],m_uiElementOrder,eIndex,paddWidth,sz,waveletC,num_wc,(double**)ws);
                            //computeCoarsenWavelets(unzippedVec[varIds[var]],offset,m_uiElementOrder,eIndex,paddWidth,sz,waveletC);
                            //l_inf=normLInfty(waveletC,NUM_COARSE_WAVELET_COEF);
                            //l_inf=normLInfty(waveletC,num_wc);
                            l_inf = normL2(waveletC,num_wc)/num_wc;
                            if(l_inf>amr_coarse_fac*tol)
                            {
                                isCoarsen=false;
                                break;
                            }

                        }


                        if(isCoarsen)
                        {

                            for(unsigned int child=0;child<NUM_CHILDREN;child++)
                            {
                                m_uiAllElements[ele+child].setFlag(((OCT_COARSE<<NUM_LEVEL_BITS)|m_uiAllElements[ele].getLevel()));
                                assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE);
                                //std::cout<<"rank: "<<m_uiActiveRank<<" element C: "<<m_uiAllElements[ele]<<" is coarsening "<<l_inf<<std::endl;
                            }

                        }


                    }

                }


            }

            delete [] waveletR;
            delete [] waveletC;
            delete [] wsIn;
            delete [] wsOut;
            delete [] ws;
            
            isOctChange=false;
            for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
                if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT) // trigger remesh only when some refinement occurs (laid back remesh :)  ) //if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)!=OCT_NO_CHANGE)
                {
                    isOctChange=true;
                    break;
                }

        }

        bool isOctChanged_g;
        MPI_Allreduce(&isOctChange,&isOctChanged_g,1,MPI_CXX_BOOL,MPI_LOR,m_uiCommGlobal);
        //if(!m_uiGlobalRank) std::cout<<"is oct changed: "<<isOctChanged_g<<std::endl;
        return isOctChanged_g;





    }

    template<typename T>
    void Mesh::getElementNodalValues(const T* vec,T* nodalValues,unsigned int elementID) const
    {

        if(!m_uiIsActive) return;
        #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                dendro::timer::t_unzip_sync_nodalval.start();
        #endif
        std::vector<T> edgeInpIn;
        std::vector<T> edgeInpOut;

        std::vector<T> faceInpIn;
        std::vector<T> faceInpOut;

        unsigned int cnum;
        bool isHanging;

        std::vector<unsigned int > edgeIndex;
        std::vector<unsigned int > faceIndex;


        bool nodeStatus[OCT_DIR_TOTAL];
        for(unsigned int w=0;w<OCT_DIR_TOTAL;w++)
            nodeStatus[w]=false;

        edgeInpIn.resize((m_uiElementOrder+1));
        edgeInpOut.resize((m_uiElementOrder+1));

        faceInpIn.resize((m_uiElementOrder+1)*(m_uiElementOrder+1));
        faceInpOut.resize((m_uiElementOrder+1)*(m_uiElementOrder+1));

        for(unsigned int k=1;k<(m_uiElementOrder);k++)
            for(unsigned int j=1;j<(m_uiElementOrder);j++)
                for(unsigned int i=1;i<(m_uiElementOrder);i++)
                {
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
                    assert(!(this->isNodeHanging(elementID,i,j,k))); // internal nodes cannot be hangging.
                }
        nodeStatus[OCT_DIR_INTERNAL]=true;

        // face interpolations
        // face : OCT_DIR_LEFT (1)
        isHanging=this->isFaceHanging(elementID,OCT_DIR_LEFT,cnum);
        if(isHanging) {
            faceNodesIndex(elementID, OCT_DIR_LEFT, faceIndex, false);
            for (unsigned int index = 0; index < faceIndex.size(); index++)
                faceInpIn[index] = vec[m_uiE2NMapping_CG[faceIndex[index]]];

            this->parent2ChildInterpolation(&(*(faceInpIn.begin())), &(*(faceInpOut.begin())), cnum, 2);

            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]=faceInpOut[k*(m_uiElementOrder+1)+j];

            nodeStatus[OCT_DIR_LEFT_DOWN]=true;
            nodeStatus[OCT_DIR_LEFT_UP]=true;
            nodeStatus[OCT_DIR_LEFT_BACK]=true;
            nodeStatus[OCT_DIR_LEFT_FRONT]=true;


            nodeStatus[OCT_DIR_LEFT_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_LEFT_UP_BACK]=true;
            nodeStatus[OCT_DIR_LEFT_UP_FRONT]=true;
            nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]=true;



        }else
        {

            for(unsigned int k=1;k<m_uiElementOrder;k++)
                for(unsigned int j=1;j<m_uiElementOrder;j++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]];
        }


        // face : OCT_DIR_RIGHT (2)
        isHanging=this->isFaceHanging(elementID,OCT_DIR_RIGHT,cnum);
        if(isHanging) {
            faceNodesIndex(elementID, OCT_DIR_RIGHT, faceIndex, false);
            for (unsigned int index = 0; index < faceIndex.size(); index++)
                faceInpIn[index] = vec[m_uiE2NMapping_CG[faceIndex[index]]];

            this->parent2ChildInterpolation(&(*(faceInpIn.begin())), &(*(faceInpOut.begin())), cnum, 2);

            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]=faceInpOut[k*(m_uiElementOrder+1)+j];

            nodeStatus[OCT_DIR_RIGHT_DOWN]=true;
            nodeStatus[OCT_DIR_RIGHT_UP]=true;
            nodeStatus[OCT_DIR_RIGHT_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_FRONT]=true;


            nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_UP_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_UP_FRONT]=true;
            nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]=true;


        }else
        {
            for(unsigned int k=1;k<m_uiElementOrder;k++)
                for(unsigned int j=1;j<m_uiElementOrder;j++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]];
        }



        // face : OCT_DIR_DOWN (3)
        isHanging=this->isFaceHanging(elementID,OCT_DIR_DOWN,cnum);
        if(isHanging) {
            faceNodesIndex(elementID, OCT_DIR_DOWN, faceIndex, false);
            for (unsigned int index = 0; index < faceIndex.size(); index++)
                faceInpIn[index] = vec[m_uiE2NMapping_CG[faceIndex[index]]];

            this->parent2ChildInterpolation(&(*(faceInpIn.begin())), &(*(faceInpOut.begin())), cnum, 2);

            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]=faceInpOut[k*(m_uiElementOrder+1)+i];


            nodeStatus[OCT_DIR_RIGHT_DOWN]=true;
            nodeStatus[OCT_DIR_LEFT_DOWN]=true;
            nodeStatus[OCT_DIR_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_DOWN_FRONT]=true;

            nodeStatus[OCT_DIR_LEFT_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]=true;
            nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]=true;


        }else
        {

            for(unsigned int k=1;k<m_uiElementOrder;k++)
                for(unsigned int i=1;i<m_uiElementOrder;i++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]];
        }


        // face : OCT_DIR_UP (4)
        isHanging=this->isFaceHanging(elementID,OCT_DIR_UP,cnum);
        if(isHanging) {
            faceNodesIndex(elementID, OCT_DIR_UP, faceIndex, false);
            for (unsigned int index = 0; index < faceIndex.size(); index++)
                faceInpIn[index] = vec[m_uiE2NMapping_CG[faceIndex[index]]];

            this->parent2ChildInterpolation(&(*(faceInpIn.begin())), &(*(faceInpOut.begin())), cnum, 2);

            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]=faceInpOut[k*(m_uiElementOrder+1)+i];

            nodeStatus[OCT_DIR_RIGHT_UP]=true;
            nodeStatus[OCT_DIR_LEFT_UP]=true;
            nodeStatus[OCT_DIR_UP_BACK]=true;
            nodeStatus[OCT_DIR_UP_FRONT]=true;


            nodeStatus[OCT_DIR_LEFT_UP_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_UP_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_UP_FRONT]=true;
            nodeStatus[OCT_DIR_LEFT_UP_FRONT]=true;

        }else
        {


            for(unsigned int k=1;k<m_uiElementOrder;k++)
                for(unsigned int i=1;i<m_uiElementOrder;i++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]];
        }



        // face : OCT_DIR_BACK (5)
        isHanging=this->isFaceHanging(elementID,OCT_DIR_BACK,cnum);
        if(isHanging) {
            faceNodesIndex(elementID, OCT_DIR_BACK, faceIndex, false);
            for (unsigned int index = 0; index < faceIndex.size(); index++)
                faceInpIn[index] = vec[m_uiE2NMapping_CG[faceIndex[index]]];

            this->parent2ChildInterpolation(&(*(faceInpIn.begin())), &(*(faceInpOut.begin())), cnum, 2);

            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=faceInpOut[j*(m_uiElementOrder+1)+i];

            nodeStatus[OCT_DIR_LEFT_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_BACK]=true;
            nodeStatus[OCT_DIR_UP_BACK]=true;
            nodeStatus[OCT_DIR_DOWN_BACK]=true;

            nodeStatus[OCT_DIR_LEFT_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_LEFT_UP_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]=true;
            nodeStatus[OCT_DIR_RIGHT_UP_BACK]=true;


        }else
        {

            for(unsigned int j=1;j<m_uiElementOrder;j++)
                for(unsigned int i=1;i<m_uiElementOrder;i++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
        }


        // face : OCT_DIR_FRONT (6)
        isHanging=this->isFaceHanging(elementID,OCT_DIR_FRONT,cnum);
        if(isHanging) {
            faceNodesIndex(elementID, OCT_DIR_FRONT, faceIndex, false);
            for (unsigned int index = 0; index < faceIndex.size(); index++)
                faceInpIn[index] = vec[m_uiE2NMapping_CG[faceIndex[index]]];

            this->parent2ChildInterpolation(&(*(faceInpIn.begin())), &(*(faceInpOut.begin())), cnum, 2);

            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=faceInpOut[j*(m_uiElementOrder+1)+i];


            nodeStatus[OCT_DIR_LEFT_FRONT]=true;
            nodeStatus[OCT_DIR_RIGHT_FRONT]=true;
            nodeStatus[OCT_DIR_UP_FRONT]=true;
            nodeStatus[OCT_DIR_DOWN_FRONT]=true;

            nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]=true;
            nodeStatus[OCT_DIR_LEFT_UP_FRONT]=true;
            nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]=true;
            nodeStatus[OCT_DIR_RIGHT_UP_FRONT]=true;


        }else
        {
            for(unsigned int j=1;j<m_uiElementOrder;j++)
                for(unsigned int i=1;i<m_uiElementOrder;i++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
        }



        // edge: OCT_DIR_LEFT_DOWN (1)

        if((!nodeStatus[OCT_DIR_LEFT_DOWN]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_LEFT_DOWN,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_LEFT,OCT_DIR_DOWN,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]=edgeInpOut[k];

                nodeStatus[OCT_DIR_LEFT_DOWN_BACK]=true;
                nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]=true;

            }else
            {
                for(unsigned int k=1;k<(m_uiElementOrder);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]];
            }
        }



        // edge: OCT_DIR_LEFT_UP (2)

        if((!nodeStatus[OCT_DIR_LEFT_UP]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_LEFT_UP,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_LEFT,OCT_DIR_UP,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]=edgeInpOut[k];

                nodeStatus[OCT_DIR_LEFT_UP_BACK]=true;
                nodeStatus[OCT_DIR_LEFT_UP_FRONT]=true;

            }else
            {
                for(unsigned int k=1;k<(m_uiElementOrder);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]];
            }
        }


        // edge: OCT_DIR_LEFT_BACK (3)

        if((!nodeStatus[OCT_DIR_LEFT_BACK]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_LEFT_BACK,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_LEFT,OCT_DIR_BACK,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]=edgeInpOut[j];

                nodeStatus[OCT_DIR_LEFT_DOWN_BACK]=true;
                nodeStatus[OCT_DIR_LEFT_UP_BACK]=true;

            }else
            {
                for(unsigned int j=1;j<(m_uiElementOrder);j++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]];
            }
        }



        // edge: OCT_DIR_LEFT_FRONT(4)

        if((!nodeStatus[OCT_DIR_LEFT_FRONT]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_LEFT_FRONT,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_LEFT,OCT_DIR_FRONT,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]=edgeInpOut[j];

                nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]=true;
                nodeStatus[OCT_DIR_LEFT_UP_FRONT]=true;

            }else
            {
                for(unsigned int j=1;j<(m_uiElementOrder);j++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+0]];
            }
        }



        // edge: OCT_DIR_RIGHT_DOWN (5)

        if((!nodeStatus[OCT_DIR_RIGHT_DOWN]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_RIGHT_DOWN,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_RIGHT,OCT_DIR_DOWN,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]=edgeInpOut[k];

                nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]=true;
                nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]=true;

            }else
            {
                for(unsigned int k=1;k<(m_uiElementOrder);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]];
            }
        }



        // edge: OCT_DIR_RIGHT_UP (6)

        if((!nodeStatus[OCT_DIR_RIGHT_UP]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_RIGHT_UP,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_RIGHT,OCT_DIR_UP,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]=edgeInpOut[k];

                nodeStatus[OCT_DIR_RIGHT_UP_BACK]=true;
                nodeStatus[OCT_DIR_RIGHT_UP_FRONT]=true;

            }else
            {
                for(unsigned int k=1;k<(m_uiElementOrder);k++)
                    nodalValues[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]];
            }
        }


        // edge: OCT_DIR_RIGHT_BACK (7)

        if((!nodeStatus[OCT_DIR_RIGHT_BACK]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_RIGHT_BACK,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_RIGHT,OCT_DIR_BACK,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]=edgeInpOut[j];

                nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]=true;
                nodeStatus[OCT_DIR_RIGHT_UP_BACK]=true;

            }else
            {
                for(unsigned int j=1;j<(m_uiElementOrder);j++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]];
            }
        }



        // edge: OCT_DIR_RIGHT_FRONT(8)

        if((!nodeStatus[OCT_DIR_RIGHT_FRONT]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_RIGHT_FRONT,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_RIGHT,OCT_DIR_FRONT,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]=edgeInpOut[j];

                nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]=true;
                nodeStatus[OCT_DIR_RIGHT_UP_FRONT]=true;

            }else
            {
                for(unsigned int j=1;j<(m_uiElementOrder);j++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+m_uiElementOrder]];
            }
        }





        // edge: OCT_DIR_DOWN_BACK (9)

        if((!nodeStatus[OCT_DIR_DOWN_BACK]))
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_DOWN_BACK,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_DOWN,OCT_DIR_BACK,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]=edgeInpOut[i];

                nodeStatus[OCT_DIR_LEFT_DOWN_BACK]=true;
                nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]=true;

            }else
            {
                for(unsigned int i=1;i<(m_uiElementOrder);i++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]];
            }
        }


        // edge: OCT_DIR_DOWN_FRONT (10)

        if(!nodeStatus[OCT_DIR_DOWN_FRONT])
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_DOWN_FRONT,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_DOWN,OCT_DIR_FRONT,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]=edgeInpOut[i];

                nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]=true;
                nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]=true;

            }else
            {
                for(unsigned int i=1;i<(m_uiElementOrder);i++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+i]];
            }
        }


        // edge: OCT_DIR_UP_BACK (11)

        if(!nodeStatus[OCT_DIR_UP_BACK])
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_UP_BACK,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_UP,OCT_DIR_BACK,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]=edgeInpOut[i];

                nodeStatus[OCT_DIR_LEFT_UP_BACK]=true;
                nodeStatus[OCT_DIR_RIGHT_UP_BACK]=true;

            }else
            {
                for(unsigned int i=1;i<(m_uiElementOrder);i++)
                    nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]];
            }
        }



        // edge: OCT_DIR_UP_FRONT (12)

        if(!nodeStatus[OCT_DIR_UP_FRONT])
        {
            isHanging=this->isEdgeHanging(elementID,OCT_DIR_UP_FRONT,cnum);
            if(isHanging)
            {
                edgeNodeIndex(elementID,OCT_DIR_UP,OCT_DIR_FRONT,edgeIndex,false);
                for(unsigned int index=0;index<edgeIndex.size();index++)
                    edgeInpIn[index]=vec[m_uiE2NMapping_CG[edgeIndex[index]]];

                this->parent2ChildInterpolation(&(*(edgeInpIn.begin())),&(*(edgeInpOut.begin())),cnum,1);

                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]=edgeInpOut[i];

                nodeStatus[OCT_DIR_LEFT_UP_FRONT]=true;
                nodeStatus[OCT_DIR_RIGHT_UP_FRONT]=true;


            }else
            {
                for(unsigned int i=1;i<(m_uiElementOrder);i++)
                    nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+i]];
            }
        }


        //node: OCT_DIR_LEFT_DOWN_BACK
        if((!(this->isNodeHanging(elementID,0,0,0))) || (!nodeStatus[OCT_DIR_LEFT_DOWN_BACK]))
            nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]];

        //node: OCT_DIR_RIGHT_DOWN_BACK
        if(!(this->isNodeHanging(elementID,m_uiElementOrder,0,0)) || (!nodeStatus[OCT_DIR_RIGHT_DOWN_BACK]))
            nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]];

        //node: OCT_DIR_LEFT_UP_BACK
        if(!(this->isNodeHanging(elementID,0,m_uiElementOrder,0)) || (!nodeStatus[OCT_DIR_LEFT_UP_BACK]))
            nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]];

        //node: OCT_DIR_RIGHT_UP_BACK
        if(!(this->isNodeHanging(elementID,m_uiElementOrder,m_uiElementOrder,0)) || (!nodeStatus[OCT_DIR_RIGHT_UP_BACK]))
            nodalValues[0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+0*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]];


        //node: OCT_DIR_LEFT_DOWN_FRONT
        if(!(this->isNodeHanging(elementID,0,0,m_uiElementOrder))|| (!nodeStatus[OCT_DIR_LEFT_DOWN_FRONT]))
            nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+0]];

        //node: OCT_DIR_RIGHT_DOWN_FRONT
        if(!(this->isNodeHanging(elementID,m_uiElementOrder,0,m_uiElementOrder))|| (!nodeStatus[OCT_DIR_RIGHT_DOWN_FRONT]))
            nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+0*(m_uiElementOrder+1)+m_uiElementOrder]];

        //node: OCT_DIR_LEFT_UP_FRONT
        if(!(this->isNodeHanging(elementID,0,m_uiElementOrder,m_uiElementOrder)) || (!nodeStatus[OCT_DIR_LEFT_UP_FRONT]))
            nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+0]];

        //node: OCT_DIR_RIGHT_UP_FRONT
        if(!(this->isNodeHanging(elementID,m_uiElementOrder,m_uiElementOrder,m_uiElementOrder)) || (!nodeStatus[OCT_DIR_RIGHT_UP_FRONT]))
            nodalValues[m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]=vec[m_uiE2NMapping_CG[elementID*m_uiNpE+m_uiElementOrder*(m_uiElementOrder+1)*(m_uiElementOrder+1)+m_uiElementOrder*(m_uiElementOrder+1)+m_uiElementOrder]];

        #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                dendro::timer::t_unzip_sync_nodalval.stop();
        #endif

    }

    template<typename T>
    void Mesh::computeElementalContribution(const T* in,T* out,unsigned int elementID) const
    {
        if(!m_uiIsActive) return;

        const unsigned int eleOrder = m_uiElementOrder;
        const unsigned int npe_1d = eleOrder + 1;
        const unsigned int npe_2d = (eleOrder + 1) * (eleOrder + 1);
        const unsigned int nPe = (eleOrder + 1) * (eleOrder + 1) * (eleOrder + 1);

        //@todo later move this to outer allocation and reuse the memeory. 
        double* qMat = new double[nPe*nPe];
        double* qTIn = new double[nPe];

        this->getElementQMat(elementID,qMat,true);
        
        for (unsigned int i=0;i<nPe; i++)
        {
            qTIn[i] = 0;

            for(unsigned int j=0;j<nPe;j++)
            {
                qTIn[i] += qMat[j*nPe + i] *in[j]; // note the transpose. 
            }
        }

        for (unsigned int i=0;i<nPe; i++)
            out[m_uiE2NMapping_CG[elementID*nPe + i]] += qTIn[i];


        delete [] qMat;
        delete [] qTIn;

        
        return;


    }

    template<typename T>
    void Mesh::interGridTransfer(std::vector<T> & vec,const ot::Mesh* pMesh)
    {


        MPI_Comm comm=m_uiCommGlobal;
        int rank,npes;

        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        int * sendNodeCount= new int[npes];
        int * recvNodeCount= new int[npes];
        int * sendNodeOffset = new int [npes];
        int * recvNodeOffset = new int [npes];
        std::vector<T> wVec; // dg of m2prime;


        for(unsigned int p=0;p<npes;p++)
            sendNodeCount[p]=0;

        if(m_uiIsActive)
        {
            MPI_Comm comm1=m_uiCommActive;
            const int rank1=m_uiActiveRank;
            const int npes1=m_uiActiveNpes;

            //1. compute the number of m2 octants (based of m1 splitters)
            unsigned int m2primeCount=0;
            for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
            {
                if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                    m2primeCount+=NUM_CHILDREN;
                else if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE)
                {
                    assert(m_uiAllElements[ele].getParent()==m_uiAllElements[ele+NUM_CHILDREN-1].getParent());
                    m2primeCount+=1;
                    ele+=(NUM_CHILDREN-1);
                }else
                {
                    assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_NO_CHANGE);
                    m2primeCount+=1;
                }

            }

            const unsigned int numM2PrimeElems=m2primeCount;

            std::vector<T> nodalVals;
            nodalVals.resize(m_uiNpE);

            std::vector<T> interp_out;
            interp_out.resize(m_uiNpE);

            wVec.resize(numM2PrimeElems*m_uiNpE);

            std::vector<ot::TreeNode> m2prime; // m2 partiioned with m1 splitters.

            m2primeCount=0;
            unsigned int cnum;
            bool isHanging;
            for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
            {
                if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                {
                    m_uiAllElements[ele].addChildren(m2prime);


                    this->getElementNodalValues(&(*(vec.begin())),&(*(nodalVals.begin())),ele);
                    for(unsigned int child=0;child<NUM_CHILDREN;child++)
                    {
                        cnum=m2prime[m2primeCount+child].getMortonIndex();
                        this->parent2ChildInterpolation(&(*(nodalVals.begin())),&(*(wVec.begin()+(m2primeCount+child)*m_uiNpE)),cnum,3);
                    }

                    m2primeCount+=NUM_CHILDREN;

                }
                else if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE)
                {
                    assert(m_uiAllElements[ele].getParent()==m_uiAllElements[ele+NUM_CHILDREN-1].getParent());
                    m2prime.push_back(m_uiAllElements[ele].getParent());

                    if(m_uiElementOrder==1)
                    {

                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                        {
                            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                                    {

                                        isHanging=this->isNodeHanging((ele+child),i,j,k);
                                        if(isHanging)
                                        {
                                            wVec[m2primeCount*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];

                                        }else{
                                            
                                            cnum=m_uiAllElements[(ele+child)].getMortonIndex();
                                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                            {
                                                wVec[ m2primeCount*m_uiNpE +  (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
                                                
                                            }

                                        }

                                    }

                        }

                    }else
                    {
                        
                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                        {
                            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                                    {

                                        isHanging=this->isNodeHanging((ele+child),i,j,k);
                                        if(isHanging)
                                        {
                                            wVec[m2primeCount*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];

                                        }else{
                                            
                                            cnum=m_uiAllElements[(ele+child)].getMortonIndex();
                                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                            {
                                                wVec[ m2primeCount*m_uiNpE +  (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
                                                
                                            }

                                        }

                                    }

                        }

                    }


                    ele+=(NUM_CHILDREN-1);
                    m2primeCount+=1;

                }else
                {
                    assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_NO_CHANGE);
                    m2prime.push_back(m_uiAllElements[ele]);

                    this->getElementNodalValues(&(*(vec.begin())),&(*(wVec.begin()+(m2primeCount*m_uiNpE))),ele);
                    m2primeCount+=1;


                }

            }


            assert(seq::test::isUniqueAndSorted(m2prime));

            if(npes1==1 && pMesh->isActive() && pMesh->getMPICommSize()==1)
            {

                // sequential case.

                if((wVec.size()/m_uiNpE)!=pMesh->getNumLocalMeshElements())
                    std::cout<<"rank1: "<<rank1<<" seq::[Inter-grid Transfer error ]: Recvn DG elements: "<<(wVec.size()/m_uiNpE)<<" m2 num local elements "<<pMesh->getNumLocalMeshElements()<<std::endl;

                assert((wVec.size()/m_uiNpE)==pMesh->getNumLocalMeshElements());

                std::vector<T>tVec;
                pMesh->createVector<T>(tVec,0);
                const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));

                const unsigned int m2LocalElemBegin=pMesh->getElementLocalBegin();
                const unsigned int m2LocalElemEnd=pMesh->getElementLocalEnd();

                const unsigned int m2LocalNodeBegin=pMesh->getNodeLocalBegin();
                const unsigned int m2LocalNodeEnd=pMesh->getNodeLocalEnd();

                unsigned int lookUp;
                const unsigned int eleOrder=pMesh->getElementOrder();

                for(unsigned int ele=m2LocalElemBegin;ele<m2LocalElemEnd;ele++)
                {
                    for(unsigned int k=0;k<eleOrder+1;k++)
                        for(unsigned int j=0;j<eleOrder+1;j++)
                            for(unsigned int i=0;i<eleOrder+1;i++)
                            {
                                if(!(pMesh->isNodeHanging(ele,i,j,k)))
                                {
                                    lookUp=e2n[ele*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                    if((lookUp>=m2LocalNodeBegin && lookUp<m2LocalNodeEnd) )
                                        tVec[lookUp]=wVec[(ele-m2LocalElemBegin)*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                }



                            }

                }


                std::swap(vec,tVec);
                tVec.clear();
                return ;
            }

            int npes2=0;
            int rank2=0;
            std::vector<ot::TreeNode> m2_splitters;
            //note : assumes that global rank 0 is going to be active always. 
            if(pMesh->isActive())
            {
                npes2=pMesh->getMPICommSize();
                rank2=pMesh->getMPIRank();
                const std::vector<ot::TreeNode> m2_splitters_root=pMesh->getSplitterElements();
                m2_splitters.resize(2*npes2);
                for(unsigned int w=0;w<m2_splitters_root.size();w++)
                    m2_splitters[w]=m2_splitters_root[w];
            } 
            
            par::Mpi_Bcast(&npes2,1,0,comm1);
            par::Mpi_Bcast(&rank2,1,0,comm1);
            m2_splitters.resize(2*npes2);
            par::Mpi_Bcast(&(*(m2_splitters.begin())),2*npes2,0,comm1);
            assert(seq::test::isUniqueAndSorted(m2_splitters));



                std::vector<ot::SearchKey> m2primeSK;
                m2primeSK.resize(m2prime.size());

                for(unsigned int e=0;e<m2prime.size();e++)
                {
                    m2primeSK[e]=ot::SearchKey(m2prime[e]);
                    m2primeSK[e].addOwner(rank1); // note that this is the rank in comm1. 
                }


                std::vector<ot::Key> m2_splitterKeys;
                m2_splitterKeys.resize(2*npes2);

                for(unsigned int p=0;p<npes2;p++)
                {
                    m2_splitterKeys[2*p]=ot::Key(m2_splitters[2*p]);
                    m2_splitterKeys[2*p].addOwner(p);

                    m2_splitterKeys[2*p+1]=ot::Key(m2_splitters[2*p+1]);
                    m2_splitterKeys[2*p+1].addOwner(p);

                    m2primeSK.push_back(ot::SearchKey(m2_splitters[2*p]));
                    m2primeSK.push_back(ot::SearchKey(m2_splitters[2*p+1]));
                }

                ot::SearchKey rootSK(m_uiDim,m_uiMaxDepth);
                std::vector<ot::SearchKey> tmpNodes;

                SFC::seqSort::SFC_treeSort(&(*(m2primeSK.begin())),m2primeSK.size(),tmpNodes,tmpNodes,tmpNodes,m_uiMaxDepth,m_uiMaxDepth,rootSK,ROOT_ROTATION,1,TS_SORT_ONLY);

                unsigned int skip=0;
                ot::SearchKey tmpSK;
                std::vector<ot::SearchKey> tmpSKVec;

                for(unsigned int e=0;e<(m2primeSK.size());e++)
                {
                    tmpSK=m2primeSK[e];
                    skip=1;
                    while(((e+skip)<m2primeSK.size()) && (m2primeSK[e]==m2primeSK[e+skip]))
                    {
                        if(m2primeSK[e+skip].getOwner()>=0){
                            tmpSK.addOwner(m2primeSK[e+skip].getOwner());
                        }
                        skip++;
                    }

                    tmpSKVec.push_back(tmpSK);
                    e+=(skip-1);

                }

                std::swap(m2primeSK,tmpSKVec);
                tmpSKVec.clear();

                assert(seq::test::isUniqueAndSorted(m2primeSK));
                assert(seq::test::isUniqueAndSorted(m2_splitterKeys));

                ot::Key rootKey(0,0,0,0,m_uiDim,m_uiMaxDepth);
                SFC::seqSearch::SFC_treeSearch(&(*(m2_splitterKeys.begin())),&(*(m2primeSK.begin())),0,m2_splitterKeys.size(),0,m2primeSK.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);



                unsigned int sBegin,sEnd,selectedRank;
                for(unsigned int p=0;p<npes2;p++)
                {
                    assert(m2_splitterKeys[2*p].getFlag() & OCT_FOUND);
                    assert(m2_splitterKeys[2*p+1].getFlag() & OCT_FOUND);

                    sBegin=m2_splitterKeys[2*p].getSearchResult();
                    sEnd=m2_splitterKeys[2*p+1].getSearchResult();
                    assert(sBegin<sEnd);
                    selectedRank=rankSelectRule(m_uiGlobalNpes,m_uiGlobalRank,npes2,p);
                    sendNodeCount[selectedRank]=sEnd-sBegin-1;

                    if(m2primeSK[sBegin].getOwner()>=0) sendNodeCount[selectedRank]++;
                    if(m2primeSK[sEnd].getOwner()>=0) sendNodeCount[selectedRank]++;

                    sendNodeCount[selectedRank]*=m_uiNpE;

                }

                // we don't need below for intergrid transfer, but these can be help full for debugging.
                m2prime.clear();
                m2primeSK.clear();




        }


        par::Mpi_Alltoall(sendNodeCount,recvNodeCount,1,comm);

        sendNodeOffset[0]=0;
        recvNodeOffset[0]=0;

        omp_par::scan(sendNodeCount,sendNodeOffset,npes);
        omp_par::scan(recvNodeCount,recvNodeOffset,npes);


        std::vector<T> wVec_m2;
        wVec_m2.resize(recvNodeOffset[npes-1]+recvNodeCount[npes-1]);

        if((wVec_m2.size()/m_uiNpE)!=pMesh->getNumLocalMeshElements())
            std::cout<<"rank: "<<rank<<" [Inter-grid Transfer error ]: Recvn DG elements: "<<(wVec_m2.size()/m_uiNpE)<<" m2 num local elements "<<pMesh->getNumLocalMeshElements()<<std::endl;

        par::Mpi_Alltoallv_sparse(&(*(wVec.begin())),sendNodeCount,sendNodeOffset,&(*(wVec_m2.begin())),recvNodeCount,recvNodeOffset,comm);

        delete [] sendNodeCount;
        delete [] recvNodeCount;
        delete [] sendNodeOffset;
        delete [] recvNodeOffset;

        std::vector<T> tVec;
        if(pMesh->isActive())
        {
            pMesh->createVector<T>(tVec,0);
            const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));

            const unsigned int m2LocalElemBegin=pMesh->getElementLocalBegin();
            const unsigned int m2LocalElemEnd=pMesh->getElementLocalEnd();

            const unsigned int m2LocalNodeBegin=pMesh->getNodeLocalBegin();
            const unsigned int m2LocalNodeEnd=pMesh->getNodeLocalEnd();

            unsigned int lookUp;
            const unsigned int eleOrder=pMesh->getElementOrder();

            for(unsigned int ele=m2LocalElemBegin;ele<m2LocalElemEnd;ele++)
            {
                for(unsigned int k=0;k<eleOrder+1;k++)
                    for(unsigned int j=0;j<eleOrder+1;j++)
                        for(unsigned int i=0;i<eleOrder+1;i++)
                        {
                            if(!(pMesh->isNodeHanging(ele,i,j,k)))
                            {
                                lookUp=e2n[ele*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                if((lookUp>=m2LocalNodeBegin && lookUp<m2LocalNodeEnd) )
                                    tVec[lookUp]=wVec_m2[(ele-m2LocalElemBegin)*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                            }



                        }

            }


            std::swap(vec,tVec);

        }

        tVec.clear();
        return ;



    }


    template<typename T>
    void Mesh::interGridTransfer(T*& vec,const ot::Mesh* pMesh)
    {


        MPI_Comm comm=m_uiCommGlobal;
        int rank,npes;

        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        int * sendNodeCount= new int[npes];
        int * recvNodeCount= new int[npes];
        int * sendNodeOffset = new int [npes];
        int * recvNodeOffset = new int [npes];
        std::vector<T> wVec; // dg of m2prime;


        for(unsigned int p=0;p<npes;p++)
            sendNodeCount[p]=0;

        if(m_uiIsActive)
        {
            MPI_Comm comm1=m_uiCommActive;
            const int rank1=m_uiActiveRank;
            const int npes1=m_uiActiveNpes;

            //1. compute the number of m2 octants (based of m1 splitters)
            unsigned int m2primeCount=0;
            for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
            {
                if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                    m2primeCount+=NUM_CHILDREN;
                else if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE)
                {
                    assert(m_uiAllElements[ele].getParent()==m_uiAllElements[ele+NUM_CHILDREN-1].getParent());
                    m2primeCount+=1;
                    ele+=(NUM_CHILDREN-1);
                }else
                {
                    assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_NO_CHANGE);
                    m2primeCount+=1;
                }

            }

            const unsigned int numM2PrimeElems=m2primeCount;

            std::vector<T> nodalVals;
            nodalVals.resize(m_uiNpE);

            std::vector<T> interp_out;
            interp_out.resize(m_uiNpE);

            wVec.resize(numM2PrimeElems*m_uiNpE);

            std::vector<ot::TreeNode> m2prime; // m2 partiioned with m1 splitters.

            m2primeCount=0;
            unsigned int cnum;
            bool isHanging;
            for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
            {
                if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                {
                    m_uiAllElements[ele].addChildren(m2prime);


                    this->getElementNodalValues(vec,&(*(nodalVals.begin())),ele);
                    for(unsigned int child=0;child<NUM_CHILDREN;child++)
                    {
                        cnum=m2prime[m2primeCount+child].getMortonIndex();
                        this->parent2ChildInterpolation(&(*(nodalVals.begin())),&(*(wVec.begin()+(m2primeCount+child)*m_uiNpE)),cnum,3);
                    }

                    m2primeCount+=NUM_CHILDREN;

                }
                else if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE)
                {
                    assert(m_uiAllElements[ele].getParent()==m_uiAllElements[ele+NUM_CHILDREN-1].getParent());
                    m2prime.push_back(m_uiAllElements[ele].getParent());


                    if(m_uiElementOrder==1)
                    {

                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                        {
                            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                                    {

                                        isHanging=this->isNodeHanging((ele+child),i,j,k);
                                        if(isHanging)
                                        {
                                            wVec[m2primeCount*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];

                                        }else{
                                            
                                            cnum=m_uiAllElements[(ele+child)].getMortonIndex();
                                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                            {
                                                wVec[ m2primeCount*m_uiNpE +  (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
                                                
                                            }

                                        }

                                    }

                        }

                    }else
                    {
                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                        {
                            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                                    {

                                        isHanging=this->isNodeHanging((ele+child),i,j,k);
                                        if(isHanging)
                                        {
                                            wVec[m2primeCount*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];

                                        }else{
                                            
                                            cnum=m_uiAllElements[(ele+child)].getMortonIndex();
                                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                            {
                                                wVec[ m2primeCount*m_uiNpE +  (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
                                                
                                            }

                                        }

                                    }

                        }

                    }

                    


                    ele+=(NUM_CHILDREN-1);
                    m2primeCount+=1;

                }else
                {
                    assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_NO_CHANGE);
                    m2prime.push_back(m_uiAllElements[ele]);

                    this->getElementNodalValues(vec,&(*(wVec.begin()+(m2primeCount*m_uiNpE))),ele);
                    m2primeCount+=1;


                }

            }


            assert(seq::test::isUniqueAndSorted(m2prime));

            if(npes1==1 && pMesh->isActive() && pMesh->getMPICommSize()==1)
            {

                // sequential case.

                if((wVec.size()/m_uiNpE)!=pMesh->getNumLocalMeshElements())
                    std::cout<<"rank1: "<<rank1<<" seq::[Inter-grid Transfer error ]: Recvn DG elements: "<<(wVec.size()/m_uiNpE)<<" m2 num local elements "<<pMesh->getNumLocalMeshElements()<<std::endl;

                assert((wVec.size()/m_uiNpE)==pMesh->getNumLocalMeshElements());

                T * tVec=pMesh->createVector<T>(0);
                const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));

                const unsigned int m2LocalElemBegin=pMesh->getElementLocalBegin();
                const unsigned int m2LocalElemEnd=pMesh->getElementLocalEnd();

                const unsigned int m2LocalNodeBegin=pMesh->getNodeLocalBegin();
                const unsigned int m2LocalNodeEnd=pMesh->getNodeLocalEnd();

                unsigned int lookUp;
                const unsigned int eleOrder=pMesh->getElementOrder();

                for(unsigned int ele=m2LocalElemBegin;ele<m2LocalElemEnd;ele++)
                {
                    for(unsigned int k=0;k<eleOrder+1;k++)
                        for(unsigned int j=0;j<eleOrder+1;j++)
                            for(unsigned int i=0;i<eleOrder+1;i++)
                            {
                                if(!(pMesh->isNodeHanging(ele,i,j,k)))
                                {
                                    lookUp=e2n[ele*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                    if((lookUp>=m2LocalNodeBegin && lookUp<m2LocalNodeEnd) )
                                        tVec[lookUp]=wVec[(ele-m2LocalElemBegin)*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                }



                            }

                }


                std::swap(vec,tVec);
                delete [] tVec;
                return ;
            }


            int npes2=0;
            int rank2=0;
            std::vector<ot::TreeNode> m2_splitters;
            //note : assumes that global rank 0 is going to be active always. 
            if(pMesh->isActive())
            {
                npes2=pMesh->getMPICommSize();
                rank2=pMesh->getMPIRank();
                const std::vector<ot::TreeNode> m2_splitters_root=pMesh->getSplitterElements();
                m2_splitters.resize(2*npes2);
                for(unsigned int w=0;w<m2_splitters_root.size();w++)
                    m2_splitters[w]=m2_splitters_root[w];
            } 
            
            par::Mpi_Bcast(&npes2,1,0,comm1);
            par::Mpi_Bcast(&rank2,1,0,comm1);
            m2_splitters.resize(2*npes2);
            par::Mpi_Bcast(&(*(m2_splitters.begin())),2*npes2,0,comm1);
            assert(seq::test::isUniqueAndSorted(m2_splitters));
           
           
               std::vector<ot::SearchKey> m2primeSK;
               m2primeSK.resize(m2prime.size());

               for(unsigned int e=0;e<m2prime.size();e++)
               {
                   m2primeSK[e]=ot::SearchKey(m2prime[e]);
                   m2primeSK[e].addOwner(rank1); // note that this is the rank in comm1. 
               }


               std::vector<ot::Key> m2_splitterKeys;
               m2_splitterKeys.resize(2*npes2);

               for(unsigned int p=0;p<npes2;p++)
               {
                   m2_splitterKeys[2*p]=ot::Key(m2_splitters[2*p]);
                   m2_splitterKeys[2*p].addOwner(p);

                   m2_splitterKeys[2*p+1]=ot::Key(m2_splitters[2*p+1]);
                   m2_splitterKeys[2*p+1].addOwner(p);

                   m2primeSK.push_back(ot::SearchKey(m2_splitters[2*p]));
                   m2primeSK.push_back(ot::SearchKey(m2_splitters[2*p+1]));
               }

               ot::SearchKey rootSK(m_uiDim,m_uiMaxDepth);
               std::vector<ot::SearchKey> tmpNodes;

               SFC::seqSort::SFC_treeSort(&(*(m2primeSK.begin())),m2primeSK.size(),tmpNodes,tmpNodes,tmpNodes,m_uiMaxDepth,m_uiMaxDepth,rootSK,ROOT_ROTATION,1,TS_SORT_ONLY);

               unsigned int skip=0;
               ot::SearchKey tmpSK;
               std::vector<ot::SearchKey> tmpSKVec;

               for(unsigned int e=0;e<(m2primeSK.size());e++)
               {
                   tmpSK=m2primeSK[e];
                   skip=1;
                   while(((e+skip)<m2primeSK.size()) && (m2primeSK[e]==m2primeSK[e+skip]))
                   {
                       if(m2primeSK[e+skip].getOwner()>=0){
                           tmpSK.addOwner(m2primeSK[e+skip].getOwner());
                       }
                       skip++;
                   }

                   tmpSKVec.push_back(tmpSK);
                   e+=(skip-1);

               }

               std::swap(m2primeSK,tmpSKVec);
               tmpSKVec.clear();

               assert(seq::test::isUniqueAndSorted(m2primeSK));
               assert(seq::test::isUniqueAndSorted(m2_splitterKeys));

               ot::Key rootKey(0,0,0,0,m_uiDim,m_uiMaxDepth);
               SFC::seqSearch::SFC_treeSearch(&(*(m2_splitterKeys.begin())),&(*(m2primeSK.begin())),0,m2_splitterKeys.size(),0,m2primeSK.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);



               unsigned int sBegin,sEnd,selectedRank;
               for(unsigned int p=0;p<npes2;p++)
               {
                   assert(m2_splitterKeys[2*p].getFlag() & OCT_FOUND);
                   assert(m2_splitterKeys[2*p+1].getFlag() & OCT_FOUND);

                   sBegin=m2_splitterKeys[2*p].getSearchResult();
                   sEnd=m2_splitterKeys[2*p+1].getSearchResult();
                   assert(sBegin<sEnd);
                   selectedRank=rankSelectRule(m_uiGlobalNpes,m_uiGlobalRank,npes2,p);
                   sendNodeCount[selectedRank]=sEnd-sBegin-1;

                   if(m2primeSK[sBegin].getOwner()>=0) sendNodeCount[selectedRank]++;
                   if(m2primeSK[sEnd].getOwner()>=0) sendNodeCount[selectedRank]++;

                   sendNodeCount[selectedRank]*=m_uiNpE;

               }

               // we don't need below for intergrid transfer, but these can be help full for debugging.
               m2prime.clear();
               m2primeSK.clear();


        }


        par::Mpi_Alltoall(sendNodeCount,recvNodeCount,1,comm);

        sendNodeOffset[0]=0;
        recvNodeOffset[0]=0;

        omp_par::scan(sendNodeCount,sendNodeOffset,npes);
        omp_par::scan(recvNodeCount,recvNodeOffset,npes);


        std::vector<T> wVec_m2;
        wVec_m2.resize(recvNodeOffset[npes-1]+recvNodeCount[npes-1]);

        if((wVec_m2.size()/m_uiNpE)!=pMesh->getNumLocalMeshElements())
            std::cout<<"rank: "<<rank<<" [Inter-grid Transfer error ]: Recvn DG elements: "<<(wVec_m2.size()/m_uiNpE)<<" m2 num local elements "<<pMesh->getNumLocalMeshElements()<<std::endl;

        par::Mpi_Alltoallv_sparse(&(*(wVec.begin())),sendNodeCount,sendNodeOffset,&(*(wVec_m2.begin())),recvNodeCount,recvNodeOffset,comm);

        delete [] sendNodeCount;
        delete [] recvNodeCount;
        delete [] sendNodeOffset;
        delete [] recvNodeOffset;

        T * tVec=NULL;
        if(pMesh->isActive())
        {
            tVec=pMesh->createVector<T>(0);
            const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));

            const unsigned int m2LocalElemBegin=pMesh->getElementLocalBegin();
            const unsigned int m2LocalElemEnd=pMesh->getElementLocalEnd();

            const unsigned int m2LocalNodeBegin=pMesh->getNodeLocalBegin();
            const unsigned int m2LocalNodeEnd=pMesh->getNodeLocalEnd();

            unsigned int lookUp;
            const unsigned int eleOrder=pMesh->getElementOrder();

            for(unsigned int ele=m2LocalElemBegin;ele<m2LocalElemEnd;ele++)
            {
                for(unsigned int k=0;k<eleOrder+1;k++)
                    for(unsigned int j=0;j<eleOrder+1;j++)
                        for(unsigned int i=0;i<eleOrder+1;i++)
                        {
                            if(!(pMesh->isNodeHanging(ele,i,j,k)))
                            {
                                lookUp=e2n[ele*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                if((lookUp>=m2LocalNodeBegin && lookUp<m2LocalNodeEnd) )
                                    tVec[lookUp]=wVec_m2[(ele-m2LocalElemBegin)*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                            }



                        }

            }


            std::swap(vec,tVec);

        }

        delete [] tVec;
        return ;







    }

    template<typename T>
    void Mesh::interGridTransferUnzip(T*& unzip, T*& vec, const ot::Mesh *pMesh)
    {

        MPI_Comm comm=m_uiCommGlobal;
        int rank,npes;

        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        int * sendNodeCount= new int[npes];
        int * recvNodeCount= new int[npes];
        int * sendNodeOffset = new int [npes];
        int * recvNodeOffset = new int [npes];
        std::vector<T> wVec; // dg of m2prime;

        for(unsigned int p=0;p<npes;p++)
            sendNodeCount[p]=0;

        
        if(m_uiIsActive)
        {

            MPI_Comm comm1=m_uiCommActive;
            const int rank1=m_uiActiveRank;
            const int npes1=m_uiActiveNpes;

            //1. compute the number of m2 octants (based of m1 splitters)
            unsigned int m2primeCount=0;
            for(unsigned int ele=m_uiElementLocalBegin;ele<m_uiElementLocalEnd;ele++)
            {
                if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                    m2primeCount+=NUM_CHILDREN;
                else if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE)
                {
                    assert(m_uiAllElements[ele].getParent()==m_uiAllElements[ele+NUM_CHILDREN-1].getParent());
                    m2primeCount+=1;
                    ele+=(NUM_CHILDREN-1);
                }else
                {
                    assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_NO_CHANGE);
                    m2primeCount+=1;
                }

            }

            const unsigned int numM2PrimeElems=m2primeCount;

            std::vector<T> nodalVals;
            nodalVals.resize(m_uiNpE);
            
            std::vector<T> interp_out;
            interp_out.resize(m_uiNpE);
            
            const unsigned int ch_1d = 2*(m_uiElementOrder) + 1;
            std::vector<T> interp_out_all;
            interp_out_all.resize(ch_1d*ch_1d*ch_1d);

            std::vector<T> nodalVals_unzip;
            

            wVec.resize(numM2PrimeElems*m_uiNpE);
            //assert(par::test::isUniqueAndSorted(m_uiLocalBlockList,comm1));

            std::vector<ot::TreeNode> m2prime; // m2 partiioned with m1 splitters.

            m2primeCount=0;
            unsigned int cnum;
            bool isHanging;
            
            for(unsigned int blk = 0 ; blk < m_uiLocalBlockList.size(); blk++ )
            {

                const ot::TreeNode blkNode=m_uiLocalBlockList[blk].getBlockNode();
                assert(blkNode.maxX()<=m_uiMeshDomain_max && blkNode.minX()>=m_uiMeshDomain_min);
                const unsigned int regLev=m_uiLocalBlockList[blk].getRegularGridLev();
                
                const unsigned int eleIndexMin = 0;
                const unsigned int eleIndexMax = (1u<<(regLev-blkNode.getLevel()))-1;
                assert(eleIndexMax>=eleIndexMin);

                const unsigned int lx=m_uiLocalBlockList[blk].getAllocationSzX();
                const unsigned int ly=m_uiLocalBlockList[blk].getAllocationSzY();
                const unsigned int lz=m_uiLocalBlockList[blk].getAllocationSzZ();
                const unsigned int offset=m_uiLocalBlockList[blk].getOffset();
                const unsigned int paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();
                const unsigned int bflag = m_uiLocalBlockList[blk].getBlkNodeFlag();

                const unsigned int u_sz[3] = { (m_uiElementOrder+1) + 2*paddWidth , (m_uiElementOrder+1) + 2*paddWidth, (m_uiElementOrder+1) + 2*paddWidth};
                nodalVals_unzip.resize(u_sz[0]*u_sz[1]*u_sz[2]);
                
                

                //std::cout<<"blk: "<<blk<<" : "<<blkNode<<"eleB: "<<m_uiLocalBlockList[blk].getLocalElementBegin()<<" eleE: "<<m_uiLocalBlockList[blk].getLocalElementEnd()<<std::endl;

                for(unsigned int ele = m_uiLocalBlockList[blk].getLocalElementBegin(); ele < m_uiLocalBlockList[blk].getLocalElementEnd(); ele++ )
                {
                    assert(m_uiAllElements[ele].getLevel()==regLev); // this is enforced by block construction

                    const unsigned int ei=(m_uiAllElements[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                    const unsigned int ej=(m_uiAllElements[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                    const unsigned int ek=(m_uiAllElements[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                    bool fdstyle =true;
                    
                    if(isBoundaryOctant(ele))
                     fdstyle=false;

                    if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                    {
                        m_uiAllElements[ele].addChildren(m2prime);
                        if(!fdstyle)
                        {
                            // do element local interpolation
                            this->getElementNodalValues(vec,nodalVals.data(),ele);
                            for(unsigned int child=0;child<NUM_CHILDREN;child++)
                            {
                                cnum=m2prime[m2primeCount+child].getMortonIndex();
                                this->parent2ChildInterpolation(&(*(nodalVals.begin())),&(*(wVec.begin()+(m2primeCount+child)*m_uiNpE)),cnum,3);
                            }

                        }else
                        {

                            this->getUnzipElementalNodalValues(unzip,blk,ele,nodalVals_unzip.data());
                            
                            // for(unsigned int k=0;k< m_uiElementOrder+1 ; k++)
                            //  for(unsigned int j=0; j< m_uiElementOrder+1; j++)
                            //   for(unsigned int i=0; i< m_uiElementOrder+1; i++)
                            //     printf("ele: %d ijk: %d,%d,%d unzip: %f  \t zip: %f\n ",ele,i,j,k,nodalVals_unzip[(k+3)*11*11+ (j+3)*11 + (i+3)],nodalVals[k*5*5+ j*5 +i]);

                            //std::cout<<"interpolation for : "<<m_uiAllElements[ele]<<std::endl;

                            m_uiRefEl.I3D_Parent2Child_FD(nodalVals_unzip.data(),interp_out_all.data(),paddWidth);
                            for(unsigned int child=0;child<NUM_CHILDREN;child++)
                            {
                                cnum=m2prime[m2primeCount+child].getMortonIndex();
                                const char bit0 = binOp::getBit(cnum, 0);
                                const char bit1 = binOp::getBit(cnum, 1);
                                const char bit2 = binOp::getBit(cnum, 2);

                                 for(unsigned int k=0; k < (m_uiElementOrder+1);  k++)
                                  for(unsigned int j=0; j < (m_uiElementOrder+1);  j++)
                                   for(unsigned int i=0; i < (m_uiElementOrder+1);  i++)
                                   {
                                       wVec[(m2primeCount+child)*m_uiNpE + k*(m_uiElementOrder+1)*(m_uiElementOrder+1) + j* (m_uiElementOrder+1) +i] = interp_out_all[(bit2*m_uiElementOrder+k)*ch_1d*ch_1d + (bit1*m_uiElementOrder+j)*ch_1d + (bit0*m_uiElementOrder+i)];
                                   }

                            }
                        
                        }

                        m2primeCount+=NUM_CHILDREN;

                    }
                    else if((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE)
                    {

                        assert(m_uiAllElements[ele].getParent()==m_uiAllElements[ele+NUM_CHILDREN-1].getParent());
                        m2prime.push_back(m_uiAllElements[ele].getParent());
                        
                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                        {
                            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                                    {

                                        isHanging=this->isNodeHanging((ele+child),i,j,k);
                                        if(isHanging)
                                        {
                                            wVec[m2primeCount*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];

                                        }else if( (i%2==0) && (j%2==0) && (k%2==0))
                                        {
                                            cnum=m_uiAllElements[(ele+child)].getMortonIndex();
                                            wVec[m2primeCount*m_uiNpE+((((cnum & 4u)>>2u)*m_uiElementOrder+k)>>1)*(m_uiElementOrder+1)*(m_uiElementOrder+1)+((((cnum & 2u)>>1u)*m_uiElementOrder+j)>>1)*(m_uiElementOrder+1)+(((((cnum & (1u)))*m_uiElementOrder+i)>>1))]=vec[m_uiE2NMapping_CG[(ele+child)*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
                                        }

                                    }

                        }
                        ele+=(NUM_CHILDREN-1);
                        m2primeCount+=1;

                    }else
                    {
                        assert((m_uiAllElements[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_NO_CHANGE);
                        m2prime.push_back(m_uiAllElements[ele]);

                        this->getElementNodalValues(vec,&(*(wVec.begin()+(m2primeCount*m_uiNpE))),ele);
                        m2primeCount+=1;
                    }

                }
            }

            assert(seq::test::isUniqueAndSorted(m2prime));

            if(npes1==1 && pMesh->isActive() && pMesh->getMPICommSize()==1)
            {

                // sequential case.

                if((wVec.size()/m_uiNpE)!=pMesh->getNumLocalMeshElements())
                    std::cout<<"rank1: "<<rank1<<" seq::[Inter-grid Transfer error ]: Recvn DG elements: "<<(wVec.size()/m_uiNpE)<<" m2 num local elements "<<pMesh->getNumLocalMeshElements()<<std::endl;

                assert((wVec.size()/m_uiNpE)==pMesh->getNumLocalMeshElements());

                T * tVec=pMesh->createVector<T>(0);
                const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));

                const unsigned int m2LocalElemBegin=pMesh->getElementLocalBegin();
                const unsigned int m2LocalElemEnd=pMesh->getElementLocalEnd();

                const unsigned int m2LocalNodeBegin=pMesh->getNodeLocalBegin();
                const unsigned int m2LocalNodeEnd=pMesh->getNodeLocalEnd();

                unsigned int lookUp;
                const unsigned int eleOrder=pMesh->getElementOrder();

                for(unsigned int ele=m2LocalElemBegin;ele<m2LocalElemEnd;ele++)
                {
                    for(unsigned int k=0;k<eleOrder+1;k++)
                        for(unsigned int j=0;j<eleOrder+1;j++)
                            for(unsigned int i=0;i<eleOrder+1;i++)
                            {
                                if(!(pMesh->isNodeHanging(ele,i,j,k)))
                                {
                                    lookUp=e2n[ele*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                    if((lookUp>=m2LocalNodeBegin && lookUp<m2LocalNodeEnd) )
                                        tVec[lookUp]=wVec[(ele-m2LocalElemBegin)*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                }



                            }

                }


                std::swap(vec,tVec);
                delete [] tVec;
                tVec==NULL;
                return ;
            }


            int npes2=0;
            int rank2=0;
            std::vector<ot::TreeNode> m2_splitters;
            //note : assumes that global rank 0 is going to be active always. 
            if(pMesh->isActive())
            {
                npes2=pMesh->getMPICommSize();
                rank2=pMesh->getMPIRank();
                const std::vector<ot::TreeNode> m2_splitters_root=pMesh->getSplitterElements();
                m2_splitters.resize(2*npes2);
                for(unsigned int w=0;w<m2_splitters_root.size();w++)
                    m2_splitters[w]=m2_splitters_root[w];
            } 
            
            par::Mpi_Bcast(&npes2,1,0,comm1);
            par::Mpi_Bcast(&rank2,1,0,comm1);
            m2_splitters.resize(2*npes2);
            par::Mpi_Bcast(&(*(m2_splitters.begin())),2*npes2,0,comm1);
            assert(seq::test::isUniqueAndSorted(m2_splitters));
           
           
               std::vector<ot::SearchKey> m2primeSK;
               m2primeSK.resize(m2prime.size());

               for(unsigned int e=0;e<m2prime.size();e++)
               {
                   m2primeSK[e]=ot::SearchKey(m2prime[e]);
                   m2primeSK[e].addOwner(rank1); // note that this is the rank in comm1. 
               }


               std::vector<ot::Key> m2_splitterKeys;
               m2_splitterKeys.resize(2*npes2);

               for(unsigned int p=0;p<npes2;p++)
               {
                   m2_splitterKeys[2*p]=ot::Key(m2_splitters[2*p]);
                   m2_splitterKeys[2*p].addOwner(p);

                   m2_splitterKeys[2*p+1]=ot::Key(m2_splitters[2*p+1]);
                   m2_splitterKeys[2*p+1].addOwner(p);

                   m2primeSK.push_back(ot::SearchKey(m2_splitters[2*p]));
                   m2primeSK.push_back(ot::SearchKey(m2_splitters[2*p+1]));
               }

               ot::SearchKey rootSK(m_uiDim,m_uiMaxDepth);
               std::vector<ot::SearchKey> tmpNodes;

               SFC::seqSort::SFC_treeSort(&(*(m2primeSK.begin())),m2primeSK.size(),tmpNodes,tmpNodes,tmpNodes,m_uiMaxDepth,m_uiMaxDepth,rootSK,ROOT_ROTATION,1,TS_SORT_ONLY);

               unsigned int skip=0;
               ot::SearchKey tmpSK;
               std::vector<ot::SearchKey> tmpSKVec;

               for(unsigned int e=0;e<(m2primeSK.size());e++)
               {
                   tmpSK=m2primeSK[e];
                   skip=1;
                   while(((e+skip)<m2primeSK.size()) && (m2primeSK[e]==m2primeSK[e+skip]))
                   {
                       if(m2primeSK[e+skip].getOwner()>=0){
                           tmpSK.addOwner(m2primeSK[e+skip].getOwner());
                       }
                       skip++;
                   }

                   tmpSKVec.push_back(tmpSK);
                   e+=(skip-1);

               }

               std::swap(m2primeSK,tmpSKVec);
               tmpSKVec.clear();

               assert(seq::test::isUniqueAndSorted(m2primeSK));
               assert(seq::test::isUniqueAndSorted(m2_splitterKeys));

               ot::Key rootKey(0,0,0,0,m_uiDim,m_uiMaxDepth);
               SFC::seqSearch::SFC_treeSearch(&(*(m2_splitterKeys.begin())),&(*(m2primeSK.begin())),0,m2_splitterKeys.size(),0,m2primeSK.size(),m_uiMaxDepth,m_uiMaxDepth,ROOT_ROTATION);



               unsigned int sBegin,sEnd,selectedRank;
               for(unsigned int p=0;p<npes2;p++)
               {
                   assert(m2_splitterKeys[2*p].getFlag() & OCT_FOUND);
                   assert(m2_splitterKeys[2*p+1].getFlag() & OCT_FOUND);

                   sBegin=m2_splitterKeys[2*p].getSearchResult();
                   sEnd=m2_splitterKeys[2*p+1].getSearchResult();
                   assert(sBegin<sEnd);
                   selectedRank=rankSelectRule(m_uiGlobalNpes,m_uiGlobalRank,npes2,p);
                   sendNodeCount[selectedRank]=sEnd-sBegin-1;

                   if(m2primeSK[sBegin].getOwner()>=0) sendNodeCount[selectedRank]++;
                   if(m2primeSK[sEnd].getOwner()>=0) sendNodeCount[selectedRank]++;

                   sendNodeCount[selectedRank]*=m_uiNpE;

               }

               // we don't need below for intergrid transfer, but these can be help full for debugging.
               m2prime.clear();
               m2primeSK.clear();


        }

        par::Mpi_Alltoall(sendNodeCount,recvNodeCount,1,comm);

        sendNodeOffset[0]=0;
        recvNodeOffset[0]=0;

        omp_par::scan(sendNodeCount,sendNodeOffset,npes);
        omp_par::scan(recvNodeCount,recvNodeOffset,npes);


        std::vector<T> wVec_m2;
        wVec_m2.resize(recvNodeOffset[npes-1]+recvNodeCount[npes-1]);

        if((wVec_m2.size()/m_uiNpE)!=pMesh->getNumLocalMeshElements())
            std::cout<<"rank: "<<rank<<" [Inter-grid Transfer error ]: Recvn DG elements: "<<(wVec_m2.size()/m_uiNpE)<<" m2 num local elements "<<pMesh->getNumLocalMeshElements()<<std::endl;

        par::Mpi_Alltoallv_sparse(&(*(wVec.begin())),sendNodeCount,sendNodeOffset,&(*(wVec_m2.begin())),recvNodeCount,recvNodeOffset,comm);

        delete [] sendNodeCount;
        delete [] recvNodeCount;
        delete [] sendNodeOffset;
        delete [] recvNodeOffset;

        T * tVec=NULL;
        if(pMesh->isActive())
        {
            tVec=pMesh->createVector<T>(0);
            const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));

            const unsigned int m2LocalElemBegin=pMesh->getElementLocalBegin();
            const unsigned int m2LocalElemEnd=pMesh->getElementLocalEnd();

            const unsigned int m2LocalNodeBegin=pMesh->getNodeLocalBegin();
            const unsigned int m2LocalNodeEnd=pMesh->getNodeLocalEnd();

            unsigned int lookUp;
            const unsigned int eleOrder=pMesh->getElementOrder();

            for(unsigned int ele=m2LocalElemBegin;ele<m2LocalElemEnd;ele++)
            {
                for(unsigned int k=0;k<eleOrder+1;k++)
                    for(unsigned int j=0;j<eleOrder+1;j++)
                        for(unsigned int i=0;i<eleOrder+1;i++)
                        {
                            if(!(pMesh->isNodeHanging(ele,i,j,k)))
                            {
                                lookUp=e2n[ele*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                                if((lookUp>=m2LocalNodeBegin && lookUp<m2LocalNodeEnd) )
                                    tVec[lookUp]=wVec_m2[(ele-m2LocalElemBegin)*m_uiNpE+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                            }



                        }

            }


            std::swap(vec,tVec);

        }

        delete [] tVec;
        tVec=NULL;
        
        wVec.clear();
        wVec_m2.clear();

    }


    template<typename T>
    void Mesh::zip(const T* unzippedVec, T* zippedVec)
    {

        if(!m_uiIsActive) return;

        ot::TreeNode blkNode;
        unsigned int ei,ej,ek;
        unsigned int regLev;
        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        unsigned int lx,ly,lz,offset,paddWidth;

        for(unsigned int blk=0;blk<m_uiLocalBlockList.size();blk++)
        {
            blkNode=m_uiLocalBlockList[blk].getBlockNode();
            regLev=m_uiLocalBlockList[blk].getRegularGridLev();

            lx=m_uiLocalBlockList[blk].getAllocationSzX();
            ly=m_uiLocalBlockList[blk].getAllocationSzY();
            lz=m_uiLocalBlockList[blk].getAllocationSzZ();
            offset=m_uiLocalBlockList[blk].getOffset();
            paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();


            for(unsigned int elem=m_uiLocalBlockList[blk].getLocalElementBegin();elem<m_uiLocalBlockList[blk].getLocalElementEnd();elem++)
            {
                ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);


                assert(pNodes[elem].getLevel()==regLev); // this is enforced by block construction


                // (1). local nodes copy. Not need to interpolate or inject values. By block construction local octants in the block has is the same level as regular grid.
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            if((m_uiE2NMapping_DG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]/m_uiNpE)==elem)
                                zippedVec[m_uiE2NMapping_CG[elem*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]]=unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)];
                        }

            }
        }

    }


    template<typename T>
    void Mesh::OCT_DIR_LEFT_DOWN_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_LEFT_DOWN-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;


        const unsigned int ei=0;
        const unsigned int ej=0;
        unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_LEFT;
        const unsigned int dir2=OCT_DIR_DOWN;
        const unsigned int dir3=OCT_DIR_FRONT;
        const unsigned int dir4=OCT_DIR_RIGHT;
        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);


        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=(paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=(m_uiElementOrder+1);
        const unsigned int jb=(m_uiElementOrder-paddWidth+1);
        const unsigned int je=(m_uiElementOrder+1);
        const unsigned int ib=(m_uiElementOrder-paddWidth+1);
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=3;
        const unsigned int cnum2=7;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;

        unsigned int bflag=blk.getBlkNodeFlag();


        //std::cout<<" lookup : "<<pNodes[lookUp]<<" blkNode: "<<blk.getBlockNode()<<std::endl;

        while(edgeCount<blkElem_1D)
        {
            ek=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {

                getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"LEFT_DOWN_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
                // note this might not be the cnum1 cnum2.
                cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()-sz,blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ek=edgeCount+1;
                if(ek<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()-sz,blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }

                edgeCount+=2;


            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"LEFT_DOWN_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];

                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }

    template<typename T>
    void Mesh::OCT_DIR_LEFT_UP_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();


        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_LEFT_UP-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;


        const unsigned int ei=0;
        const unsigned int ej=blkElem_1D-1;
        unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_LEFT;
        const unsigned int dir2=OCT_DIR_UP;
        const unsigned int dir3=OCT_DIR_FRONT;
        const unsigned int dir4=OCT_DIR_RIGHT;
        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);


        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=(paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=(m_uiElementOrder+1);
        const unsigned int jb=0;
        const unsigned int je=paddWidth;
        const unsigned int ib=(m_uiElementOrder-paddWidth+1);
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=1;
        const unsigned int cnum2=5;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;


        while(edgeCount<blkElem_1D)
        {
            ek=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"LEFT_UP_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
                cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.maxY(),blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ek=edgeCount+1;

                if(ek<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.maxY(),blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }

                }

                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"LEFT_UP_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }



                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];

                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_LEFT_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {
        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_LEFT_BACK-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int ei=0;
        unsigned int ej=0;
        const unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_LEFT;
        const unsigned int dir2=OCT_DIR_BACK;
        const unsigned int dir3=OCT_DIR_UP;
        const unsigned int dir4=OCT_DIR_RIGHT;
        const unsigned int dir5=OCT_DIR_FRONT;
        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);

        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);




        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=paddWidth;
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth+1);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=0;
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=(m_uiElementOrder-paddWidth+1);
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=5;
        const unsigned int cnum2=7;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;


        while(edgeCount<blkElem_1D)
        {
            ej=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"LEFT_BACK_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
                cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()+ej*sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ej=edgeCount+1;
                if(ej<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()+ej*sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;


            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"LEFT_FRONT_BACK_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }

    template<typename T>
    void Mesh::OCT_DIR_LEFT_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_LEFT_FRONT-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int ei=0;
        unsigned int ej=0;
        const unsigned int ek=(blkElem_1D-1);

        const unsigned int dir1=OCT_DIR_LEFT;
        const unsigned int dir2=OCT_DIR_FRONT;
        const unsigned int dir3=OCT_DIR_UP;
        const unsigned int dir4=OCT_DIR_RIGHT;

        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);

        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=paddWidth;
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=(paddWidth);

        const unsigned int jb=0;
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=(m_uiElementOrder-paddWidth+1);
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=1;
        const unsigned int cnum2=3;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ej=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"LEFT_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
                cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()+ej*sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ej=edgeCount+1;

                if(ej<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()+ej*sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }



                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"LEFT_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }
                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_DOWN_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.



        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_RIGHT_DOWN-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int ei=(blkElem_1D-1);
        const unsigned int ej=0;
        unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_RIGHT;
        const unsigned int dir2=OCT_DIR_DOWN;
        const unsigned int dir3=OCT_DIR_FRONT;
        const unsigned int dir4=OCT_DIR_UP;
        ot::TreeNode blkNode=blk.getBlockNode();
        unsigned int sz=1u<<(m_uiMaxDepth-regLev);


        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=(paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=(m_uiElementOrder+1);
        const unsigned int jb=(m_uiElementOrder-paddWidth+1);
        const unsigned int je=(m_uiElementOrder+1);
        const unsigned int ib=0;
        const unsigned int ie=paddWidth;

        const unsigned int cnum1=2;
        const unsigned int cnum2=6;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;

        unsigned int bflag=blk.getBlkNodeFlag();


        //std::cout<<" lookup : "<<pNodes[lookUp]<<" blkNode: "<<blk.getBlockNode()<<std::endl;

        while(edgeCount<blkElem_1D)
        {
            ek=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {

                getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"RIGHT_DOWN_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
                cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()-sz,blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ek=edgeCount+1;
                if(ek<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()-sz,blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"RIGHT_DOWN_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_UP_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_RIGHT_UP-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;


        const unsigned int ei=blkElem_1D-1;
        const unsigned int ej=blkElem_1D-1;
        unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_RIGHT;
        const unsigned int dir2=OCT_DIR_UP;
        const unsigned int dir3=OCT_DIR_FRONT;
        const unsigned int dir4=OCT_DIR_UP;
        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);


        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=(paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=(m_uiElementOrder+1);
        const unsigned int jb=0;
        const unsigned int je=paddWidth;
        const unsigned int ib=0;
        const unsigned int ie=paddWidth;

        const unsigned int cnum1=0;
        const unsigned int cnum2=4;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ek=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"RIGHT_UP_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
                cnum=ot::TreeNode(blkNode.maxX(),blkNode.maxY(),blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ek=edgeCount+1;
                if(ek<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.maxX(),blkNode.maxY(),blkNode.minZ()+ek*sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }

                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"RIGHT_UP_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }

    template<typename T>
    void Mesh::OCT_DIR_RIGHT_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {
        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_RIGHT_BACK-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int ei=blkElem_1D-1;
        unsigned int ej=0;
        const unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_RIGHT;
        const unsigned int dir2=OCT_DIR_BACK;
        const unsigned int dir3=OCT_DIR_UP;
        const unsigned int dir4=OCT_DIR_FRONT;
        const unsigned int dir5=OCT_DIR_FRONT;
        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);

        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=paddWidth;
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth+1);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=0;
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=0;
        const unsigned int ie=paddWidth;

        const unsigned int cnum1=4;
        const unsigned int cnum2=6;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ej=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"LEFT_BACK_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
                cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()+ej*sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ej=edgeCount+1;
                if(ej<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()+ej*sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"LEFT_FRONT_BACK_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_RIGHT_FRONT-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int ei=(blkElem_1D-1);
        unsigned int ej=0;
        const unsigned int ek=(blkElem_1D-1);

        const unsigned int dir1=OCT_DIR_RIGHT;
        const unsigned int dir2=OCT_DIR_FRONT;
        const unsigned int dir3=OCT_DIR_UP;
        const unsigned int dir4=OCT_DIR_RIGHT;

        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);


        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=paddWidth;
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=(paddWidth);

        const unsigned int jb=0;
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=0;
        const unsigned int ie=paddWidth;

        const unsigned int cnum1=0;
        const unsigned int cnum2=2;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ej=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"LEFT_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
                cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()+ej*sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ej=edgeCount+1;
                if(ej<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()+ej*sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"LEFT_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                           else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }

    template<typename T>
    void Mesh::OCT_DIR_DOWN_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_DOWN_BACK-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;


        unsigned int ei=0;
        const unsigned int ej=0;
        const unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_DOWN;
        const unsigned int dir2=OCT_DIR_BACK;
        const unsigned int dir3=OCT_DIR_RIGHT;
        const unsigned int dir4=OCT_DIR_LEFT;
        const unsigned int dir5=OCT_DIR_UP;

        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);


        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=paddWidth;
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth+1);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=(m_uiElementOrder-paddWidth+1);
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=0;
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=6;
        const unsigned int cnum2=7;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ei=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"DOWN_BACK_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
                cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.minY()-sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ei=edgeCount+1;
                if(ei<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.minY()-sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"DOWN_BACK_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_DOWN_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_DOWN_FRONT-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;


        unsigned int ei=0;
        const unsigned int ej=0;
        const unsigned int ek=(blkElem_1D-1);

        const unsigned int dir1=OCT_DIR_DOWN;
        const unsigned int dir2=OCT_DIR_FRONT;
        const unsigned int dir3=OCT_DIR_RIGHT;
        const unsigned int dir4=OCT_DIR_UP;

        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);

        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=paddWidth;
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=paddWidth;

        const unsigned int jb=(m_uiElementOrder-paddWidth+1);
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=0;
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=2;
        const unsigned int cnum2=3;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ei=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"DOWN_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
                cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.minY()-sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ei=edgeCount+1;
                if(ei<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.minY()-sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"DOWN_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_UP_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.

        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_UP_BACK-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;


        unsigned int ei=0;
        const unsigned int ej=(blkElem_1D-1);
        const unsigned int ek=0;

        const unsigned int dir1=OCT_DIR_UP;
        const unsigned int dir2=OCT_DIR_BACK;
        const unsigned int dir3=OCT_DIR_RIGHT;
        const unsigned int dir4=OCT_DIR_LEFT;
        const unsigned int dir5=OCT_DIR_UP;

        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);

        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=paddWidth;
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth+1);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=0;
        const unsigned int je=paddWidth;

        const unsigned int ib=0;
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=4;
        const unsigned int cnum2=5;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ei=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"UP_BACK_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
                cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.maxY(),blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ei=edgeCount+1;
                if(ei<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.maxY(),blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }
                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"UP_BACK_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_UP_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {

        const unsigned int * blk2diagMap=blk.getBlk2DiagMap();
        unsigned int lookUp; // first OCT_DIR_LEFT_DOWN element.


        const unsigned int rank=getMPIRank();
        const unsigned int regLev=blk.getRegularGridLev();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());
        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int offset=blk.getOffset();
        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        unsigned int edgeCount=0;
        const unsigned int edgeDir=(OCT_DIR_UP_FRONT-EDGE_OFFSET);
        lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        unsigned int ei=0;
        const unsigned int ej=(blkElem_1D-1);
        const unsigned int ek=(blkElem_1D-1);

        const unsigned int dir1=OCT_DIR_UP;
        const unsigned int dir2=OCT_DIR_FRONT;
        const unsigned int dir3=OCT_DIR_RIGHT;
        const unsigned int dir4=OCT_DIR_UP;

        const ot::TreeNode blkNode=blk.getBlockNode();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);

        const ot::TreeNode * pNodes=&(*(m_uiAllElements.begin()));

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        const  int i_offset=paddWidth;
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=paddWidth;

        const unsigned int jb=0;
        const unsigned int je=paddWidth;

        const unsigned int ib=0;
        const unsigned int ie=(m_uiElementOrder+1);

        const unsigned int cnum1=0;
        const unsigned int cnum2=1;
        unsigned int cnum;

        unsigned int nodeLookUp_CG;
        bool isHanging;



        while(edgeCount<blkElem_1D)
        {
            ei=edgeCount;
            lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount];
            assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
            if(pNodes[lookUp].getLevel()==regLev)
            {
                getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;


            }else if(pNodes[lookUp].getLevel()<regLev)
            {
                if((pNodes[lookUp].getLevel()+1)!=regLev)
                    std::cout<<"UP_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
                cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.maxY(),blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                ei=edgeCount+1;
                if(ei<blkElem_1D)
                {
                    cnum=ot::TreeNode(blkNode.minX()+ei*sz,blkNode.maxY(),blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();
                    parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
                    for(unsigned int k=kb;k<ke;k++)
                        for(unsigned int j=jb;j<je;j++)
                            for(unsigned int i=ib;i<ie;i++)
                            {
                                unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            }

                }


                edgeCount+=2;




            }else
            {
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                if(pNodes[lookUp].getLevel()!=(regLev+1))
                    std::cout<<"DOWN_UP_FRONT_DIAG_UNIZIP ERROR: 2:1 balance error "<<std::endl;

                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum1);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }


                        }

                lookUp=blk2diagMap[edgeDir*(2*blkElem_1D)+2*edgeCount+1];
                assert(lookUp!=LOOK_UP_TABLE_DEFAULT);
                assert(pNodes[lookUp].getLevel()==(regLev+1));
                cnum=pNodes[lookUp].getMortonIndex();
                assert(cnum==cnum2);
                for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                    for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                        for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                        {
                            isHanging=isNodeHanging(lookUp,i,j,k);
                            nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            if(isHanging)
                            {
                                interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                            }
                            else
                            {
                                const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                                const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                                //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                                if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                {
                                    interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                    
                                }

                            }

                        }

                for(unsigned int k=kb;k<ke;k++)
                    for(unsigned int j=jb;j<je;j++)
                        for(unsigned int i=ib;i<ie;i++)
                        {
                            unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        }

                edgeCount+=1;



            }
        }


    }


    template<typename T>
    void Mesh::OCT_DIR_LEFT_DOWN_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_LEFT_DOWN_BACK;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=0;
        const unsigned int ej=0;
        const unsigned int ek=0;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth+1);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=(m_uiElementOrder-paddWidth+1);
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=(m_uiElementOrder-paddWidth+1);
        const unsigned int ie=(m_uiElementOrder+1);


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [LEFT_DOWN_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
            cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()-sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [LEFT_DOWN_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }
                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_DOWN_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_RIGHT_DOWN_BACK;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=blkElem_1D-1;
        const unsigned int ej=0;
        const unsigned int ek=0;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth+1);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=(m_uiElementOrder-paddWidth+1);
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=0;
        const unsigned int ie=paddWidth;


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [RIGHT_DOWN_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
            cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()-sz,blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [RIGHT_DOWN_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }

    template<typename T>
    void Mesh::OCT_DIR_LEFT_UP_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_LEFT_UP_BACK;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=0;
        const unsigned int ej=blkElem_1D-1;
        const unsigned int ek=0;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth+1);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=0;
        const unsigned int je=paddWidth;

        const unsigned int ib=(m_uiElementOrder-paddWidth+1);
        const unsigned int ie=(m_uiElementOrder+1);


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [LEFT_UP_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
            cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.maxY(),blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [LEFT_UP_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_UP_BACK_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_RIGHT_UP_BACK;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=blkElem_1D-1;
        const unsigned int ej=blkElem_1D-1;
        const unsigned int ek=0;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=-(m_uiElementOrder-paddWidth);

        const unsigned int kb=(m_uiElementOrder-paddWidth+1);
        const unsigned int ke=(m_uiElementOrder+1);

        const unsigned int jb=0;
        const unsigned int je=paddWidth;

        const unsigned int ib=0;
        const unsigned int ie=paddWidth;


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [RIGHT_UP_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
            cnum=ot::TreeNode(blkNode.maxX(),blkNode.maxY(),blkNode.minZ()-sz,regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [RIGHT_UP_BACK Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }


    template<typename T>
    void Mesh::OCT_DIR_LEFT_DOWN_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_LEFT_DOWN_FRONT;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int ei=0;
        const unsigned int ej=0;
        const unsigned int ek=blkElem_1D-1;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=paddWidth;

        const unsigned int jb=(m_uiElementOrder-paddWidth+1);
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=(m_uiElementOrder-paddWidth+1);
        const unsigned int ie=(m_uiElementOrder+1);


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [LEFT_DOWN_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
            cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.minY()-sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [LEFT_DOWN_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_DOWN_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_RIGHT_DOWN_FRONT;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=blkElem_1D-1;
        const unsigned int ej=0;
        const unsigned int ek=blkElem_1D-1;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=-(m_uiElementOrder-paddWidth);
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=paddWidth;

        const unsigned int jb=(m_uiElementOrder-paddWidth+1);
        const unsigned int je=(m_uiElementOrder+1);

        const unsigned int ib=0;
        const unsigned int ie=paddWidth;


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [RIGHT_DOWN_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
            cnum=ot::TreeNode(blkNode.maxX(),blkNode.minY()-sz,blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [RIGHT_DOWN_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }

    template<typename T>
    void Mesh::OCT_DIR_LEFT_UP_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_LEFT_UP_FRONT;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=0;
        const unsigned int ej=blkElem_1D-1;
        const unsigned int ek=blkElem_1D-1;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=-(m_uiElementOrder-paddWidth);
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=paddWidth;

        const unsigned int jb=0;
        const unsigned int je=paddWidth;

        const unsigned int ib=(m_uiElementOrder-paddWidth+1);
        const unsigned int ie=(m_uiElementOrder+1);


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [LEFT_UP_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
            cnum=ot::TreeNode(blkNode.minX()-sz,blkNode.maxY(),blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [LEFT_UP_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }


    template<typename T>
    void Mesh::OCT_DIR_RIGHT_UP_FRONT_Unzip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {
        const unsigned int rank=getMPIRank();
        const unsigned int dir=OCT_DIR_RIGHT_UP_FRONT;
        const unsigned int * blk2VertexMap=blk.getBlk2VertexMap();
        const unsigned int lookUp=blk2VertexMap[dir-VERTEX_OFFSET];

        if(lookUp==LOOK_UP_TABLE_DEFAULT) return;

        const unsigned int regLev=blk.getRegularGridLev();
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        const unsigned int offset=blk.getOffset();
        const unsigned int blkElem_1D=1u<<(regLev-blk.getBlockNode().getLevel());

        const unsigned int paddWidth=blk.get1DPadWidth();
        const unsigned int sz=1u<<(m_uiMaxDepth-regLev);
        const ot::TreeNode blkNode=blk.getBlockNode();

        const unsigned int ei=blkElem_1D-1;
        const unsigned int ej=blkElem_1D-1;
        const unsigned int ek=blkElem_1D-1;

        const unsigned int lx=blk.getAllocationSzX();
        const unsigned int ly=blk.getAllocationSzY();
        const unsigned int lz=blk.getAllocationSzZ();

        const  int i_offset=(m_uiElementOrder+paddWidth);
        const  int j_offset=(m_uiElementOrder+paddWidth);
        const  int k_offset=(m_uiElementOrder+paddWidth);

        const unsigned int kb=0;
        const unsigned int ke=paddWidth;

        const unsigned int jb=0;
        const unsigned int je=paddWidth;

        const unsigned int ib=0;
        const unsigned int ie=paddWidth;


        unsigned int cnum;
        bool isHanging;
        unsigned int nodeLookUp_CG;

        std::vector<T> interpIn;
        interpIn.resize(m_uiNpE);

        std::vector<T> interpOut;
        interpOut.resize(m_uiNpE);


        if(pNodes[lookUp].getLevel()==regLev)
        {
            getElementNodalValues(zippedVec,&(*(interpOut.begin())),lookUp);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }

        }else if(pNodes[lookUp].getLevel()<regLev)
        {
            if(pNodes[lookUp].getLevel()!=(regLev-1)) {std::cout<<"rank: "<<rank<<" [RIGHT_UP_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert(pNodes[lookUp].getLevel()==(regLev-1));
            getElementNodalValues(zippedVec,&(*(interpIn.begin())),lookUp);
            cnum=ot::TreeNode(blkNode.maxX(),blkNode.maxY(),blkNode.maxZ(),regLev,m_uiDim,m_uiMaxDepth).getMortonIndex();

            parent2ChildInterpolation(&(*(interpIn.begin())),&(*(interpOut.begin())),cnum,3);
            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }


        }else
        {
            if((pNodes[lookUp].getLevel())!=regLev+1) {std::cout<<"rank: "<<rank<<" [RIGHT_UP_FRONT Unzip]: 2:1 balance violation blk node: "<<blkNode<<" lookup : "<<pNodes[lookUp]<<std::endl; exit(0);}
            assert((pNodes[lookUp].getLevel())==regLev+1);

            cnum=pNodes[lookUp].getMortonIndex();
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
                for(unsigned int j=0;j<m_uiElementOrder+1;j++)
                    for(unsigned int i=0;i<m_uiElementOrder+1;i++)
                    {
                        isHanging=isNodeHanging(lookUp,i,j,k);
                        nodeLookUp_CG=m_uiE2NMapping_CG[lookUp*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                        if(isHanging)
                        {
                            interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=zippedVec[nodeLookUp_CG];
                        }
                        else
                        {
                            const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                            const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                            const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                            //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                            if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                            {
                                interpOut[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = zippedVec[nodeLookUp_CG];
                                
                            }

                        }

                    }


            for(unsigned int k=kb;k<ke;k++)
                for(unsigned int j=jb;j<je;j++)
                    for(unsigned int i=ib;i<ie;i++)
                    {
                        unzippedVec[offset+(ek*m_uiElementOrder+k+k_offset)*(ly*lx)+(ej*m_uiElementOrder+j+j_offset)*(lx)+(ei*m_uiElementOrder+i+i_offset)]=interpOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                    }



        }




    }


    template<typename T>
    void Mesh::blockDiagonalUnZip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {

        // NOTE!!!!!: When copying the diagnal padding skip the first and last point of the block padding. Since they contains the advective derivterm.
        // This is not the correct value when at the case 3 of diagonal padding in all edge directions.
        OCT_DIR_LEFT_DOWN_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_LEFT_UP_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_LEFT_BACK_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_LEFT_FRONT_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_RIGHT_DOWN_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_RIGHT_UP_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_RIGHT_BACK_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_RIGHT_FRONT_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_DOWN_BACK_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_DOWN_FRONT_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_UP_BACK_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_UP_FRONT_Unzip(blk,zippedVec,unzippedVec);


    }

    template<typename T>
    void Mesh::blockVertexUnZip(const ot::Block & blk,const T* zippedVec, T* unzippedVec)
    {

        OCT_DIR_LEFT_DOWN_BACK_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_RIGHT_DOWN_BACK_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_LEFT_UP_BACK_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_RIGHT_UP_BACK_Unzip(blk,zippedVec,unzippedVec);

        OCT_DIR_LEFT_DOWN_FRONT_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_RIGHT_DOWN_FRONT_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_LEFT_UP_FRONT_Unzip(blk,zippedVec,unzippedVec);
        OCT_DIR_RIGHT_UP_FRONT_Unzip(blk,zippedVec,unzippedVec);
    }

    template<typename T>
    void Mesh::child2ParentInjection(const T* in , T* out, unsigned int* child, unsigned int lev) const 
    {

        for(unsigned int cnum=0;cnum<NUM_CHILDREN;cnum++)
        {
            if(child[cnum] == LOOK_UP_TABLE_DEFAULT || m_uiAllElements[child[cnum]].getLevel()!= lev ||  !m_uiIsNodalMapValid[child[cnum]]) continue;
            
            for(unsigned int k=0;k<m_uiElementOrder+1;k++)
            for(unsigned int j=0;j<m_uiElementOrder+1;j++)
            for(unsigned int i=0;i<m_uiElementOrder+1;i++)
            {
                const bool isHanging=this->isNodeHanging(child[cnum],i,j,k);
                if(isHanging)
                {
                    out[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]=in[m_uiE2NMapping_CG[child[cnum]*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i]];
                }
                else
                {
                    const unsigned int iix = m_uiElementOrder * (int) (cnum & 1u)  +  i;
                    const unsigned int jjy = m_uiElementOrder * (int) ((cnum & 2u)>>1u)  +  j;
                    const unsigned int kkz = m_uiElementOrder * (int) ((cnum & 4u)>>2u)  +  k;
                    //std::cout<<" iix: "<<iix<<" jjy: "<<jjy<<" kkz: "<<kkz<<std::endl;

                    if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                    {
                        out[ (kkz>>1u) * (m_uiElementOrder+1)*(m_uiElementOrder+1) + (jjy>>1u) * (m_uiElementOrder+1)+(iix>>1u)] = in [m_uiE2NMapping_CG [child[cnum]*m_uiNpE+k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j *(m_uiElementOrder+1)+ i ] ];
                        
                    }

                }

            }
        }

    }


    template <typename T>
    void Mesh::unzip(const T* zippedVec, T* unzippedVec)
    {
        if(!m_uiIsActive) return;

        if(m_uiLocalBlockList.empty())
            return;


        ot::TreeNode blkNode;
        unsigned int ei,ej,ek; // element wise xyz coordinates.
        const ot::TreeNode * pNodes= &(*(m_uiAllElements.begin()));
        unsigned int regLev;
        //unsigned int blkNpe_1D;

        unsigned int lookUp;
        unsigned int lookUp1;
        unsigned int cnum;
        unsigned int faceCnum;

        unsigned int faceNeighCnum1[4]={0,0,0,0}; // immidiate neighbors
        unsigned int faceNeighCnum2[4]={0,0,0,0}; // neighbor's neighbors


        register unsigned int nodeLookUp_CG;
        register unsigned int nodeLookUp_DG;

        std::vector<T> interpOrInjectionOut; // interpolation or injection output.
        std::vector<T> injectionInput;// input for the injection (values from all the 8 children) (This should be put in the order of the morton ordering. )
        std::vector<T> interpolationInput;

        std::vector<T> edgeInterpIn;
        std::vector<T> edgeInterpOut;

        std::vector<T> faceInterpIn;
        std::vector<T> faceInterpOut;

        std::vector<T> lookUpElementVec;
        std::vector<T> parentEleInterpIn;
        std::vector<T> parentEleInterpOut;

        std::vector<unsigned int> edgeIndex;
        std::vector<unsigned int > faceIndex;
        std::vector<unsigned int > child;
        child.resize(NUM_CHILDREN);

        interpOrInjectionOut.resize(m_uiNpE);
        interpolationInput.resize(m_uiNpE);
        //injectionInput.resize(m_uiNpE*NUM_CHILDREN);

        std::vector<T> injectionTest;
        injectionTest.resize(m_uiNpE*NUM_CHILDREN);

        edgeIndex.resize((m_uiElementOrder+1));
        faceIndex.resize((m_uiElementOrder+1)*(m_uiElementOrder+1));

        edgeInterpIn.resize((m_uiElementOrder+1));
        edgeInterpOut.resize((m_uiElementOrder+1));

        faceInterpIn.resize((m_uiElementOrder+1)*(m_uiElementOrder+1));
        faceInterpOut.resize((m_uiElementOrder+1)*(m_uiElementOrder+1));

        lookUpElementVec.resize(m_uiNpE);
        parentEleInterpIn.resize(m_uiNpE);
        parentEleInterpOut.resize(m_uiNpE);

        unsigned int mid_bit=0;
        unsigned int sz;
        bool isHanging;
        unsigned int ownerID,ii_x,jj_y,kk_z;
        unsigned int eleIndexMin=0;
        unsigned int eleIndexMax=0;
        bool edgeHanging;
        bool faceHanging;

        unsigned int lx,ly,lz,offset,paddWidth;
        bool isParentValue=false;

        unsigned int fid[(NUM_CHILDREN>>1u)];
        unsigned int cid[(NUM_CHILDREN>>1u)];

        /*if(!rank) std::cout<<"begin unzip "<<std::endl;*/
        #ifdef DEBUG_UNZIP_OP
            double d_min,d_max;
            d_min=-0.5;
            d_max=0.5;
            double x,y,z;
            unsigned int x1,y1,z1;
            std::function<double(double,double,double)> func =[d_min,d_max](const double x,const double y,const double z){ return (sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
        #endif

        // NOTE: Be careful when you access ghost elements for padding. (You should only access the level 1 ghost elements. You should not access the level 2 ghost elements at any time. )
        paddWidth = m_uiLocalBlockList[0].get1DPadWidth();
        if(m_uiElementOrder ==4 && paddWidth==3)
        {
            //std::cout<<"read spt points : "<<m_uiElementOrder<<" pwidth : "<<paddWidth<<std::endl;
            readSpecialPtsBegin(zippedVec);
        }
            

        for(unsigned int blk=0;blk<m_uiLocalBlockList.size();blk++)
        {
            blkNode=m_uiLocalBlockList[blk].getBlockNode();
            assert(blkNode.maxX()<=m_uiMeshDomain_max && blkNode.minX()>=m_uiMeshDomain_min);
            regLev=m_uiLocalBlockList[blk].getRegularGridLev();
            //blkNpe_1D=m_uiElementOrder*(1u<<(regLev-blkNode.getLevel()))+1+2*GHOST_WIDTH;
            //std::cout<<"rank: "<<m_uiActiveRank<<" -- blkNpw_1D: "<<blkNpe_1D<<" blkNode: "<<blkNode<<" regLev: "<<regLev<<std::endl;

            sz=1u<<(m_uiMaxDepth-regLev);
            eleIndexMax=(1u<<(regLev-blkNode.getLevel()))-1;
            assert(eleIndexMax>=eleIndexMin);

            lx=m_uiLocalBlockList[blk].getAllocationSzX();
            ly=m_uiLocalBlockList[blk].getAllocationSzY();
            lz=m_uiLocalBlockList[blk].getAllocationSzZ();
            offset=m_uiLocalBlockList[blk].getOffset();
            paddWidth=m_uiLocalBlockList[blk].get1DPadWidth();


            for(unsigned int elem=m_uiLocalBlockList[blk].getLocalElementBegin();elem<m_uiLocalBlockList[blk].getLocalElementEnd();elem++)
            {
                ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                //std::cout<<"blk: "<<blk<<" : "<<blkNode<<" ek: "<<(ek)<<" ej: "<<(ej)<<" ei: "<<(ei)<<" elem: "<<m_uiAllElements[elem]<<std::endl;
                assert(pNodes[elem].getLevel()==regLev); // this is enforced by block construction
                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_internal.start();
                #endif
                this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),elem);
                //this->getElementNodalValues(zippedVec,&(*(parentEleInterpIn.begin())),elem);
                // note: do not change the parentInterpIn values. These are used to interpolate the 3rd point in the advective terms.
                for(unsigned int w=0;w<m_uiNpE;w++)
                    parentEleInterpIn[w]=lookUpElementVec[w];

                // (1). local nodes copy. Not need to interpolate or inject values. By block construction local octants in the block has is the same level as regular grid.
                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_cpy.start();
                #endif
                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                        for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                           unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_cpy.stop();
                    dendro::timer::t_unzip_sync_internal.stop();
                #endif
                // (2). copy the ghost layer (we only copy GHOST_WIDTH amounts of data from the zipped array )z`
                //---------------------------------------------------------X direction padding --------------------------------------------------------------------------------------------------------------------
                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                            dendro::timer::t_unzip_sync_face[0].start();
                #endif
                if((pNodes[elem].minX()==blkNode.minX()))
                {
                    assert(ei==eleIndexMin);

                    lookUp=m_uiE2EMapping[elem*m_uiNumDirections+OCT_DIR_LEFT];
                    if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                    {

                        if(pNodes[lookUp].getLevel()==pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c1.start();
                            #endif
                            assert(paddWidth<(m_uiElementOrder+1));
                            this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),lookUp);

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif
                            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                    for(unsigned int i=(m_uiElementOrder-paddWidth);i<(m_uiElementOrder+1);i++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i-(m_uiElementOrder-paddWidth))]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c1.stop();
                            #endif

                        }else if(pNodes[lookUp].getLevel()<pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c2.start();
                            #endif
                            assert(pNodes[lookUp].getLevel()+1==regLev);
                            mid_bit=m_uiMaxDepth - pNodes[lookUp].getLevel()-1;
                            cnum=( ((((pNodes[elem].getZ()) >> mid_bit) & 1u) << 2u) | ((((pNodes[elem].getY()) >> mid_bit) & 1u) << 1u) | ((((pNodes[elem].getX()-sz)) >>mid_bit) & 1u));
                            //std::cout<<"elem: "<<elem<<" : "<<m_uiAllElements[elem]<<" lookup: "<<m_uiAllElements[lookUp]<<" child: "<<ot::TreeNode(pNodes[elem].getX()-sz,pNodes[elem].getY(),pNodes[elem].getZ(),pNodes[elem].getLevel(),m_uiDim,m_uiMaxDepth)<<" cnum: "<<cnum<<std::endl;


                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                 dendro::timer::t_unzip_sync_cpy.start();
                             #endif


                            #ifdef USE_FD_INTERP_FOR_UNZIP
                                const int st = this->getBlkBdyParentCNums(blk,elem,OCT_DIR_LEFT,child.data(),fid,cid);
                                if(st > 0)
                                {
                                    const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
                                    this->getBlkBoundaryParentNodes(zippedVec, lookUpElementVec.data(), interpolationInput.data(), interpOrInjectionOut.data(), lookUp, fid, cid,child.data());
                                    for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
                                    {
                                        assert(pNodes[lookUp] == pNodes[m_uiE2EMapping[child[fid[w]]*m_uiNumDirections + OCT_DIR_LEFT]]);
                                        assert(child[fid[w]] != LOOK_UP_TABLE_DEFAULT);
                                        //assert(m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin) ] == blk);

                                        if(child[fid[w]]<m_uiElementLocalBegin || child[fid[w]]>=m_uiElementLocalEnd)
                                            continue;

                                        this->parent2ChildInterpolation(lookUpElementVec.data(),interpOrInjectionOut.data(),cid[w],m_uiDim);


                                        const ot::Block blk_fd = m_uiLocalBlockList[m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin)]];
                                        const ot::TreeNode blkNode_fd = blk_fd.getBlockNode();
                                        const unsigned int regL_fd = blk_fd.getRegularGridLev();
                                        
                                        const unsigned int lx_fd = blk_fd.getAllocationSzX();
                                        const unsigned int ly_fd = blk_fd.getAllocationSzY();
                                        const unsigned int lz_fd = blk_fd.getAllocationSzZ();

                                        const unsigned int offset_fd = blk_fd.getOffset();


                                        
                                        const unsigned int ei_fd = (pNodes[child[fid[w]]].getX()-blkNode_fd.getX())>>(m_uiMaxDepth-regL_fd);
                                        const unsigned int ej_fd = (pNodes[child[fid[w]]].getY()-blkNode_fd.getY())>>(m_uiMaxDepth-regL_fd);
                                        const unsigned int ek_fd = (pNodes[child[fid[w]]].getZ()-blkNode_fd.getZ())>>(m_uiMaxDepth-regL_fd);

                                        assert(paddWidth<(m_uiElementOrder+1));
                                        for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                        for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                        for(unsigned int i=(m_uiElementOrder-paddWidth);i<(m_uiElementOrder+1);i++)
                                            unzippedVec[offset_fd+(ek_fd*m_uiElementOrder+k+paddWidth)*(ly_fd*lx_fd)+(ej_fd*m_uiElementOrder+j+paddWidth)*(lx_fd)+(ei_fd*m_uiElementOrder+i-(m_uiElementOrder-paddWidth))]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                    }

                                    


                                }
                            #else
                                this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),lookUp);
                                this->parent2ChildInterpolation(&(*(lookUpElementVec.begin())),&(*(interpOrInjectionOut.begin())),cnum);

                                
                                assert(paddWidth<(m_uiElementOrder+1));
                                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                        for(unsigned int i=(m_uiElementOrder-paddWidth);i<(m_uiElementOrder+1);i++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i-(m_uiElementOrder-paddWidth))]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            #endif

                            

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop() ;
                                dendro::timer::t_unzip_sync_f_c2.stop();
                            #endif


                        }else if(pNodes[lookUp].getLevel()>pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c3.start();
                            #endif
                            assert(pNodes[lookUp].getLevel()==(regLev+1));
                            //child.resize(NUM_CHILDREN,LOOK_UP_TABLE_DEFAULT);
                            // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                            child[1]=lookUp;
                            child[3]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_UP];
                            assert(child[3]!=LOOK_UP_TABLE_DEFAULT);
                            child[5]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_FRONT];
                            assert(child[5]!=LOOK_UP_TABLE_DEFAULT);
                            child[7]=m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_FRONT];
                            assert(child[7]!=LOOK_UP_TABLE_DEFAULT);

                            if(m_uiElementOrder ==4 && paddWidth==3)
                            {
                                // we need to search for the additional points. 
                                child[0]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_LEFT];
                                child[2]=m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_LEFT];
                                child[4]=m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_LEFT];
                                child[6]=m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_LEFT];

                            }else
                            {
                                child[0]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_LEFT];
                                child[2]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_LEFT];
                                child[4]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_LEFT];
                                child[6]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_LEFT];

                            }


                            this->child2ParentInjection(zippedVec,interpOrInjectionOut.data(),child.data(),pNodes[lookUp].getLevel());

                            

                            

                            #ifdef DEBUG_UNZIP_OP_3PT
                                faceNeighCnum1[0]=1;faceNeighCnum1[1]=3;faceNeighCnum1[2]=5;faceNeighCnum1[3]=7;
                                faceNeighCnum2[0]=0;faceNeighCnum2[1]=2;faceNeighCnum2[2]=4;faceNeighCnum2[3]=6;



                                for(unsigned int index=0;index<(NUM_CHILDREN>>1u);index++)
                                {
                                    interpDownWind(fd::D1_ORDER_4_DOWNWIND,elem,child[faceNeighCnum1[index]],&(*(lookUpElementVec.begin())),faceNeighCnum2[index],&(*(parentEleInterpIn.begin())),&(*(parentEleInterpOut.begin())),OCT_DIR_LEFT,paddWidth,zippedVec,&(*(interpOrInjectionOut.begin())));
                                }
                            #endif

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif
                            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                    for(unsigned int i=(m_uiElementOrder-paddWidth);i<(m_uiElementOrder+1);i++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i-(m_uiElementOrder-paddWidth))]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c3.stop();
                            #endif


                        }

                    }

                }

                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_face[0].stop();
                #endif

                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_face[1].start();
                #endif

                if((pNodes[elem].maxX()==blkNode.maxX()))
                {
                    assert(ei==eleIndexMax);
                    lookUp=m_uiE2EMapping[elem*m_uiNumDirections+OCT_DIR_RIGHT];
                    if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                    {
                        if(pNodes[lookUp].getLevel()==pNodes[elem].getLevel())
                        {
                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c1.start();
                            #endif
                            assert(paddWidth<(m_uiElementOrder+1));
                            this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),lookUp);

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif
                            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                    for(unsigned int i=0;i<(paddWidth+1);i++)
                                      unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+((ei+1)*m_uiElementOrder+paddWidth+i)]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c1.stop();
                            #endif


                        }else if(pNodes[lookUp].getLevel()<pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c2.start();
                            #endif
                            assert(pNodes[lookUp].getLevel()+1==regLev);
                            mid_bit=m_uiMaxDepth - pNodes[lookUp].getLevel()-1;
                            cnum=( ((((pNodes[elem].getZ()) >> mid_bit) & 1u) << 2u) | ((((pNodes[elem].getY()) >> mid_bit) & 1u) << 1u) | ((((pNodes[elem].getX()+sz)) >>mid_bit) & 1u));
                            //std::cout<<"elem: "<<elem<<" : "<<m_uiAllElements[elem]<<" lookup: "<<m_uiAllElements[lookUp]<<" child: "<<ot::TreeNode(pNodes[elem].getX()+sz,pNodes[elem].getY(),pNodes[elem].getZ(),pNodes[elem].getLevel(),m_uiDim,m_uiMaxDepth)<<" cnum: "<<cnum<<std::endl;

                            

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif

                            #ifdef USE_FD_INTERP_FOR_UNZIP
                                const int st = this->getBlkBdyParentCNums(blk,elem,OCT_DIR_RIGHT,child.data(),fid,cid);
                                if(st > 0)
                                {
                                    const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
                                    this->getBlkBoundaryParentNodes(zippedVec, lookUpElementVec.data(), interpolationInput.data(), interpOrInjectionOut.data(), lookUp, fid, cid,child.data());
                                    for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
                                    {
                                        assert(pNodes[lookUp] == pNodes[m_uiE2EMapping[child[fid[w]]*m_uiNumDirections + OCT_DIR_RIGHT]]);
                                        assert(child[fid[w]] != LOOK_UP_TABLE_DEFAULT);

                                        if(child[fid[w]]<m_uiElementLocalBegin || child[fid[w]]>=m_uiElementLocalEnd)
                                            continue;

                                        this->parent2ChildInterpolation(lookUpElementVec.data(),interpOrInjectionOut.data(),cid[w],m_uiDim);
                                        
                                        const ot::Block blk_fd = m_uiLocalBlockList[m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin)]];
                                        const ot::TreeNode blkNode_fd = blk_fd.getBlockNode();
                                        const unsigned int regL_fd = blk_fd.getRegularGridLev();
                                        
                                        const unsigned int lx_fd = blk_fd.getAllocationSzX();
                                        const unsigned int ly_fd = blk_fd.getAllocationSzY();
                                        const unsigned int lz_fd = blk_fd.getAllocationSzZ();

                                        const unsigned int offset_fd = blk_fd.getOffset();


                                        
                                        const unsigned int ei_fd = (pNodes[child[fid[w]]].getX()-blkNode_fd.getX())>>(m_uiMaxDepth-regL_fd);
                                        const unsigned int ej_fd = (pNodes[child[fid[w]]].getY()-blkNode_fd.getY())>>(m_uiMaxDepth-regL_fd);
                                        const unsigned int ek_fd = (pNodes[child[fid[w]]].getZ()-blkNode_fd.getZ())>>(m_uiMaxDepth-regL_fd);


                                        for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                        for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                            for(unsigned int i=0;i<(paddWidth+1);i++)
                                            unzippedVec[offset_fd+(ek_fd*m_uiElementOrder+k+paddWidth)*(ly_fd*lx_fd)+(ej_fd*m_uiElementOrder+j+paddWidth)*(lx_fd)+((ei_fd+1)*m_uiElementOrder+paddWidth+i)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                    }

                                    


                                }
                            #else
                                this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),lookUp);
                                this->parent2ChildInterpolation(&(*(lookUpElementVec.begin())),&(*(interpOrInjectionOut.begin())),cnum);
                                assert(paddWidth<(m_uiElementOrder+1));
                                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                        for(unsigned int i=0;i<(paddWidth+1);i++)
                                            unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+((ei+1)*m_uiElementOrder+paddWidth+i)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            #endif



                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c2.stop();
                            #endif


                        }else if(pNodes[lookUp].getLevel()>pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c3.start();
                            #endif
                            // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                            child[0]=lookUp;
                            child[2]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_UP];
                            assert(child[2]!=LOOK_UP_TABLE_DEFAULT);
                            child[4]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_FRONT];
                            assert(child[4]!=LOOK_UP_TABLE_DEFAULT);
                            child[6]=m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_FRONT];
                            assert(child[6]!=LOOK_UP_TABLE_DEFAULT);

                            if(m_uiElementOrder ==4 && paddWidth==3)
                            {
                                child[1]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_RIGHT];
                                child[3]=m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_RIGHT];
                                child[5]=m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_RIGHT];
                                child[7]=m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_RIGHT];

                            }else
                            {
                                child[1]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_RIGHT];
                                child[3]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_RIGHT];
                                child[5]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_RIGHT];
                                child[7]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_RIGHT];

                            }

                            


                            this->child2ParentInjection(zippedVec,interpOrInjectionOut.data(),child.data(),pNodes[lookUp].getLevel());

                            #ifdef DEBUG_UNZIP_OP_3PT
                                faceNeighCnum1[0]=0;faceNeighCnum1[1]=2;faceNeighCnum1[2]=4;faceNeighCnum1[3]=6;
                                faceNeighCnum2[0]=1;faceNeighCnum2[1]=3;faceNeighCnum2[2]=5;faceNeighCnum2[3]=7;


                                for(unsigned int index=0;index<(NUM_CHILDREN>>1u);index++)
                                {
                                    interpUpWind(fd::D1_ORDER_4_UPWIND,elem,child[faceNeighCnum1[index]],&(*(lookUpElementVec.begin())),faceNeighCnum2[index],&(*(parentEleInterpIn.begin())),&(*(parentEleInterpOut.begin())),OCT_DIR_RIGHT,paddWidth,zippedVec,&(*(interpOrInjectionOut.begin())));
                                }
                            #endif

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif
                            
                            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                    for(unsigned int i=0;i<(paddWidth+1);i++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+((ei+1)*m_uiElementOrder+paddWidth+i)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c3.stop();
                            #endif


                        }

                    }

                }

                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_face[1].stop();
                #endif
                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_face[2].start();
                #endif

                //--------------------------------------------------------------------------------------------------- Y Direction----------------------------------------------------------------------------------
                if((pNodes[elem].minY()==blkNode.minY()))
                {
                    assert(ej==0);

                    lookUp=m_uiE2EMapping[elem*m_uiNumDirections+OCT_DIR_DOWN];
                    if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                    {

                        if(pNodes[lookUp].getLevel()==pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c1.start();
                            #endif
                            assert(paddWidth<(m_uiElementOrder+1));
                            this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),lookUp);

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif
                            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int j=(m_uiElementOrder-paddWidth);j<(m_uiElementOrder+1);j++)
                                       unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j-(m_uiElementOrder-paddWidth))*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c1.stop();
                            #endif

                        }else if(pNodes[lookUp].getLevel()<pNodes[elem].getLevel())
                        {
                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c2.start();
                            #endif
                            assert(pNodes[lookUp].getLevel()+1==regLev);
                            mid_bit=m_uiMaxDepth - pNodes[lookUp].getLevel()-1;
                            cnum=( ((((pNodes[elem].getZ()) >> mid_bit) & 1u) << 2u) | (((((pNodes[elem].getY()-sz)) >> mid_bit) & 1u) << 1u) | (((pNodes[elem].getX()) >>mid_bit) & 1u));

                            //std::cout<<"elem: "<<elem<<" : "<<m_uiAllElements[elem]<<" lookup: "<<m_uiAllElements[lookUp]<<" child: "<<ot::TreeNode(pNodes[elem].getX()-sz,pNodes[elem].getY(),pNodes[elem].getZ(),pNodes[elem].getLevel(),m_uiDim,m_uiMaxDepth)<<" cnum: "<<cnum<<std::endl;
                            

                            //std::cout<<"m_uiActiveRank : "<<m_uiActiveRank<<"parent to child interpolation executed"<<std::endl;
                            assert(paddWidth<(m_uiElementOrder+1));
                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif

                            #ifdef USE_FD_INTERP_FOR_UNZIP
                            const int st = this->getBlkBdyParentCNums(blk,elem,OCT_DIR_DOWN,child.data(),fid,cid);
                            if(st > 0)
                            {
                                const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
                                this->getBlkBoundaryParentNodes(zippedVec, lookUpElementVec.data(), interpolationInput.data(), interpOrInjectionOut.data(), lookUp, fid, cid,child.data());
                                for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
                                {
                                    assert(pNodes[lookUp] == pNodes[m_uiE2EMapping[child[fid[w]]*m_uiNumDirections + OCT_DIR_DOWN]]);
                                    assert(child[fid[w]] != LOOK_UP_TABLE_DEFAULT);

                                    if(child[fid[w]]<m_uiElementLocalBegin || child[fid[w]]>=m_uiElementLocalEnd)
                                        continue;
                                        
                                    this->parent2ChildInterpolation(lookUpElementVec.data(),interpOrInjectionOut.data(),cid[w],m_uiDim);

                                    const ot::Block blk_fd = m_uiLocalBlockList[m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin)]];
                                    const ot::TreeNode blkNode_fd = blk_fd.getBlockNode();
                                    const unsigned int regL_fd = blk_fd.getRegularGridLev();
                                    
                                    const unsigned int lx_fd = blk_fd.getAllocationSzX();
                                    const unsigned int ly_fd = blk_fd.getAllocationSzY();
                                    const unsigned int lz_fd = blk_fd.getAllocationSzZ();

                                    const unsigned int offset_fd = blk_fd.getOffset();

                                    const unsigned int ei_fd = (pNodes[child[fid[w]]].getX()-blkNode_fd.getX())>>(m_uiMaxDepth-regL_fd);
                                    const unsigned int ej_fd = (pNodes[child[fid[w]]].getY()-blkNode_fd.getY())>>(m_uiMaxDepth-regL_fd);
                                    const unsigned int ek_fd = (pNodes[child[fid[w]]].getZ()-blkNode_fd.getZ())>>(m_uiMaxDepth-regL_fd);
                                    
                                    for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                        for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                            for(unsigned int j=(m_uiElementOrder-paddWidth);j<(m_uiElementOrder+1);j++)
                                                unzippedVec[offset_fd+(ek_fd*m_uiElementOrder+k+paddWidth)*(ly_fd*lx_fd)+(ej_fd*m_uiElementOrder+j-(m_uiElementOrder-paddWidth))*(lx_fd)+(ei_fd*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                }

                                


                            }
                            #else
                                this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),lookUp);
                                this->parent2ChildInterpolation(&(*(lookUpElementVec.begin())),&(*(interpOrInjectionOut.begin())),cnum);

                                for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int j=(m_uiElementOrder-paddWidth);j<(m_uiElementOrder+1);j++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j-(m_uiElementOrder-paddWidth))*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];


                            #endif


                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c2.stop();
                            #endif


                        }else if(pNodes[lookUp].getLevel()>pNodes[elem].getLevel())
                        {
                            // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c3.start();
                            #endif
                            child[2]=lookUp;
                            child[3]=m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_RIGHT];
                            assert(child[3]!=LOOK_UP_TABLE_DEFAULT);
                            child[6]=m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_FRONT];
                            assert(child[6]!=LOOK_UP_TABLE_DEFAULT);
                            child[7]=m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_FRONT];
                            assert(child[7]!=LOOK_UP_TABLE_DEFAULT);

                            if(m_uiElementOrder ==4 && paddWidth==3)
                            {
                                child[0]=m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_DOWN];
                                child[1]=m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_DOWN];
                                child[4]=m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_DOWN];
                                child[5]=m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_DOWN];


                            }else
                            {
                                child[0]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_DOWN];
                                child[1]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_DOWN];
                                child[4]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_DOWN];
                                child[5]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_DOWN];

                            }

                            
                            this->child2ParentInjection(zippedVec,interpOrInjectionOut.data(),child.data(),pNodes[lookUp].getLevel());

                            #ifdef DEBUG_UNZIP_OP_3PT
                                faceNeighCnum1[0]=2;faceNeighCnum1[1]=3;faceNeighCnum1[2]=6;faceNeighCnum1[3]=7;
                                faceNeighCnum2[0]=0;faceNeighCnum2[1]=1;faceNeighCnum2[2]=4;faceNeighCnum2[3]=5;


                                for(unsigned int index=0;index<(NUM_CHILDREN>>1u);index++)
                                {
                                interpDownWind(fd::D1_ORDER_4_DOWNWIND,elem,child[faceNeighCnum1[index]],&(*(lookUpElementVec.begin())),faceNeighCnum2[index],&(*(parentEleInterpIn.begin())),&(*(parentEleInterpOut.begin())),OCT_DIR_DOWN,paddWidth,zippedVec,&(*(interpOrInjectionOut.begin())));
                                }
                            #endif


                            //std::cout<<"m_uiActiveRank : "<<m_uiActiveRank<<"child to parent interpolation executed"<<std::endl;
                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif
                            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int j=(m_uiElementOrder-paddWidth);j<(m_uiElementOrder+1);j++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+(ej*m_uiElementOrder+j-(m_uiElementOrder-paddWidth))*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c3.stop();
                            #endif

                        }

                    }

                }

                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_face[2].stop();
                #endif
                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_face[3].start();
                #endif
                if((pNodes[elem].maxY()==blkNode.maxY()))
                {
                    assert(ej==(1u<<(regLev-blkNode.getLevel()))-1);
                    lookUp=m_uiE2EMapping[elem*m_uiNumDirections+OCT_DIR_UP];
                    if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                    {
                        if(pNodes[lookUp].getLevel()==pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c1.start();
                            #endif
                            assert(paddWidth<(m_uiElementOrder+1));
                            this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),lookUp);

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif
                            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int j=0;j<(paddWidth+1);j++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+((ej+1)*m_uiElementOrder+paddWidth+j)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c1.stop();
                            #endif


                        }else if(pNodes[lookUp].getLevel()<pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c2.start();
                            #endif
                            assert(pNodes[lookUp].getLevel()+1==regLev);
                            mid_bit=m_uiMaxDepth - pNodes[lookUp].getLevel()-1;
                            cnum=( ((((pNodes[elem].getZ()) >> mid_bit) & 1u) << 2u) | (((((pNodes[elem].getY()+sz)) >> mid_bit) & 1u) << 1u) | (((pNodes[elem].getX()) >>mid_bit) & 1u));
                            //std::cout<<"elem: "<<elem<<" : "<<m_uiAllElements[elem]<<" lookup: "<<m_uiAllElements[lookUp]<<" child: "<<ot::TreeNode(pNodes[elem].getX()+sz,pNodes[elem].getY(),pNodes[elem].getZ(),pNodes[elem].getLevel(),m_uiDim,m_uiMaxDepth)<<" cnum: "<<cnum<<std::endl;
                            

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif

                            

                            #ifdef USE_FD_INTERP_FOR_UNZIP
                            const int st = this->getBlkBdyParentCNums(blk,elem,OCT_DIR_UP,child.data(),fid,cid);
                            if(st > 0)
                            {
                                const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
                                this->getBlkBoundaryParentNodes(zippedVec, lookUpElementVec.data(), interpolationInput.data(), interpOrInjectionOut.data(), lookUp, fid, cid,child.data());
                                for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
                                {
                                    assert(pNodes[lookUp] == pNodes[m_uiE2EMapping[child[fid[w]]*m_uiNumDirections + OCT_DIR_UP]]);
                                    assert(child[fid[w]] != LOOK_UP_TABLE_DEFAULT);

                                    if(child[fid[w]]<m_uiElementLocalBegin || child[fid[w]]>=m_uiElementLocalEnd)
                                        continue;

                                    this->parent2ChildInterpolation(lookUpElementVec.data(),interpOrInjectionOut.data(),cid[w],m_uiDim);


                                    const ot::Block blk_fd = m_uiLocalBlockList[m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin)]];
                                    const ot::TreeNode blkNode_fd = blk_fd.getBlockNode();
                                    const unsigned int regL_fd = blk_fd.getRegularGridLev();
                                    
                                    const unsigned int lx_fd = blk_fd.getAllocationSzX();
                                    const unsigned int ly_fd = blk_fd.getAllocationSzY();
                                    const unsigned int lz_fd = blk_fd.getAllocationSzZ();

                                    const unsigned int offset_fd = blk_fd.getOffset();
                                    
                                    const unsigned int ei_fd = (pNodes[child[fid[w]]].getX()-blkNode_fd.getX())>>(m_uiMaxDepth-regL_fd);
                                    const unsigned int ej_fd = (pNodes[child[fid[w]]].getY()-blkNode_fd.getY())>>(m_uiMaxDepth-regL_fd);
                                    const unsigned int ek_fd = (pNodes[child[fid[w]]].getZ()-blkNode_fd.getZ())>>(m_uiMaxDepth-regL_fd);
                                    
                                    for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                        for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                            for(unsigned int j=0;j<(paddWidth+1);j++)
                                                unzippedVec[offset_fd+(ek_fd*m_uiElementOrder+k+paddWidth)*(ly_fd*lx_fd)+((ej_fd+1)*m_uiElementOrder+paddWidth+j)*(lx_fd)+(ei_fd*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                }

                                


                            }
                            #else
                            this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),lookUp);
                            this->parent2ChildInterpolation(&(*(lookUpElementVec.begin())),&(*(interpOrInjectionOut.begin())),cnum);

                            assert(paddWidth<(m_uiElementOrder+1));
                            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int j=0;j<(paddWidth+1);j++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+((ej+1)*m_uiElementOrder+paddWidth+j)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #endif


                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c2.stop();
                            #endif



                        }else if(pNodes[lookUp].getLevel()>pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c3.start();
                            #endif
                            // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                            child[0]=lookUp;
                            child[1]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_RIGHT];
                            assert(child[1]!=LOOK_UP_TABLE_DEFAULT);
                            child[4]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_FRONT];
                            assert(child[4]!=LOOK_UP_TABLE_DEFAULT);
                            child[5]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_FRONT];
                            assert(child[5]!=LOOK_UP_TABLE_DEFAULT);

                            if(m_uiElementOrder ==4 && paddWidth==3)
                            {
                                child[2]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_UP];
                                child[3]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_UP];
                                child[6]=m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_UP];
                                child[7]=m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_UP];


                            }else
                            {
                                child[2]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_UP];
                                child[3]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_UP];
                                child[6]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_UP];
                                child[7]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_UP];
                            }

                            


                            this->child2ParentInjection(zippedVec,interpOrInjectionOut.data(),child.data(),pNodes[lookUp].getLevel());

                            #ifdef DEBUG_UNZIP_OP_3PT
                                faceNeighCnum1[0]=0;faceNeighCnum1[1]=1;faceNeighCnum1[2]=4;faceNeighCnum1[3]=5;
                                faceNeighCnum2[0]=2;faceNeighCnum2[1]=3;faceNeighCnum2[2]=6;faceNeighCnum2[3]=7;


                                for(unsigned int index=0;index<(NUM_CHILDREN>>1u);index++)
                                {
                                interpUpWind(fd::D1_ORDER_4_UPWIND,elem,child[faceNeighCnum1[index]],&(*(lookUpElementVec.begin())),faceNeighCnum2[index],&(*(parentEleInterpIn.begin())),&(*(parentEleInterpOut.begin())),OCT_DIR_UP,paddWidth,zippedVec,&(*(interpOrInjectionOut.begin())));
                                }
                            #endif

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif
                            for(unsigned int k=0;k<(m_uiElementOrder+1);k++)
                                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int j=0;j<(paddWidth+1);j++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k+paddWidth)*(ly*lx)+((ej+1)*m_uiElementOrder+paddWidth+j)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];


                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c3.stop();
                            #endif

                        }

                    }

                }

                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_face[3].stop();
                #endif
                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_face[4].start();
                #endif
                //--------------------------------------------------------------------- Z direction padding. -------------------------------------------------------------------------------------------------------

                if((pNodes[elem].minZ()==blkNode.minZ()))
                {
                    assert(ek==0);

                    lookUp=m_uiE2EMapping[elem*m_uiNumDirections+OCT_DIR_BACK];
                    if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                    {

                        if(pNodes[lookUp].getLevel()==pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c1.start();
                            #endif
                            assert(paddWidth<(m_uiElementOrder+1));
                            this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),lookUp);

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif
                            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int k=(m_uiElementOrder-paddWidth);k<(m_uiElementOrder+1);k++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k-(m_uiElementOrder-paddWidth))*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c1.stop();
                            #endif


                        }else if(pNodes[lookUp].getLevel()<pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c2.start();
                            #endif
                            assert(pNodes[lookUp].getLevel()+1==regLev);
                            mid_bit=m_uiMaxDepth - pNodes[lookUp].getLevel()-1;
                            cnum=( (((((pNodes[elem].getZ()-sz)) >> mid_bit) & 1u) << 2u) | ((((pNodes[elem].getY()) >> mid_bit) & 1u) << 1u) | (((pNodes[elem].getX()) >>mid_bit) & 1u));
                            //std::cout<<"elem: "<<elem<<" : "<<m_uiAllElements[elem]<<" lookup: "<<m_uiAllElements[lookUp]<<" child: "<<ot::TreeNode(pNodes[elem].getX()-sz,pNodes[elem].getY(),pNodes[elem].getZ(),pNodes[elem].getLevel(),m_uiDim,m_uiMaxDepth)<<" cnum: "<<cnum<<std::endl;
                            
                            //std::cout<<"m_uiActiveRank : "<<m_uiActiveRank<<"parent to child interpolation executed"<<std::endl;
                            assert(paddWidth<(m_uiElementOrder+1));
                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif

                            

                            #ifdef USE_FD_INTERP_FOR_UNZIP
                            const int st = this->getBlkBdyParentCNums(blk,elem,OCT_DIR_BACK,child.data(),fid,cid);
                            if(st > 0)
                            {
                                const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
                                this->getBlkBoundaryParentNodes(zippedVec, lookUpElementVec.data(), interpolationInput.data(), interpOrInjectionOut.data(), lookUp, fid, cid,child.data());
                                for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
                                {
                                    assert(pNodes[lookUp] == pNodes[m_uiE2EMapping[child[fid[w]]*m_uiNumDirections + OCT_DIR_BACK]]);
                                    assert(child[fid[w]] != LOOK_UP_TABLE_DEFAULT);

                                    if(child[fid[w]]<m_uiElementLocalBegin || child[fid[w]]>=m_uiElementLocalEnd)
                                        continue;
                                    
                                    this->parent2ChildInterpolation(lookUpElementVec.data(),interpOrInjectionOut.data(),cid[w],m_uiDim);
                                    const ot::Block blk_fd = m_uiLocalBlockList[m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin)]];
                                    const ot::TreeNode blkNode_fd = blk_fd.getBlockNode();
                                    const unsigned int regL_fd = blk_fd.getRegularGridLev();
                                    
                                    const unsigned int lx_fd = blk_fd.getAllocationSzX();
                                    const unsigned int ly_fd = blk_fd.getAllocationSzY();
                                    const unsigned int lz_fd = blk_fd.getAllocationSzZ();

                                    const unsigned int offset_fd = blk_fd.getOffset();


                                    
                                    const unsigned int ei_fd = (pNodes[child[fid[w]]].getX()-blkNode_fd.getX())>>(m_uiMaxDepth-regL_fd);
                                    const unsigned int ej_fd = (pNodes[child[fid[w]]].getY()-blkNode_fd.getY())>>(m_uiMaxDepth-regL_fd);
                                    const unsigned int ek_fd = (pNodes[child[fid[w]]].getZ()-blkNode_fd.getZ())>>(m_uiMaxDepth-regL_fd);
                                    

                                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int k=(m_uiElementOrder-paddWidth);k<(m_uiElementOrder+1);k++)
                                        unzippedVec[offset_fd+(ek_fd*m_uiElementOrder+k-(m_uiElementOrder-paddWidth))*(ly_fd*lx_fd)+(ej_fd*m_uiElementOrder+j+paddWidth)*(lx_fd)+(ei_fd*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                                }

                                


                            }
                            #else
                            this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),lookUp);
                            this->parent2ChildInterpolation(&(*(lookUpElementVec.begin())),&(*(interpOrInjectionOut.begin())),cnum);

                            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int k=(m_uiElementOrder-paddWidth);k<(m_uiElementOrder+1);k++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k-(m_uiElementOrder-paddWidth))*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #endif
                            

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c2.stop();
                            #endif


                        }else if(pNodes[lookUp].getLevel()>pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c3.start();
                            #endif
                            // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                            child[4]=lookUp;
                            child[5]=m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_RIGHT];
                            assert(child[5]!=LOOK_UP_TABLE_DEFAULT);
                            child[6]=m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_UP];
                            assert(child[6]!=LOOK_UP_TABLE_DEFAULT);
                            child[7]=m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_UP];
                            assert(child[7]!=LOOK_UP_TABLE_DEFAULT);

                            if(m_uiElementOrder ==4 && paddWidth==3)
                            {
                                child[0]=m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_BACK];
                                child[1]=m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_BACK];
                                child[2]=m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_BACK];
                                child[3]=m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_BACK];

                            }else
                            {
                                child[0]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[4]*m_uiNumDirections+OCT_DIR_BACK];
                                child[1]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[5]*m_uiNumDirections+OCT_DIR_BACK];
                                child[2]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[6]*m_uiNumDirections+OCT_DIR_BACK];
                                child[3]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[7]*m_uiNumDirections+OCT_DIR_BACK];
                            }

                            

                            this->child2ParentInjection(zippedVec,interpOrInjectionOut.data(),child.data(),pNodes[lookUp].getLevel());


                            #ifdef DEBUG_UNZIP_OP_3PT
                                faceNeighCnum1[0]=4;faceNeighCnum1[1]=5;faceNeighCnum1[2]=6;faceNeighCnum1[3]=7;
                                faceNeighCnum2[0]=0;faceNeighCnum2[1]=1;faceNeighCnum2[2]=2;faceNeighCnum2[3]=3;


                                for(unsigned int index=0;index<(NUM_CHILDREN>>1u);index++)
                                {
                                    interpDownWind(fd::D1_ORDER_4_DOWNWIND,elem,child[faceNeighCnum1[index]],&(*(lookUpElementVec.begin())),faceNeighCnum2[index],&(*(parentEleInterpIn.begin())),&(*(parentEleInterpOut.begin())),OCT_DIR_BACK,paddWidth,zippedVec,&(*(interpOrInjectionOut.begin())));
                                }
                            #endif



                            //std::cout<<"m_uiActiveRank : "<<m_uiActiveRank<<"child to parent interpolation executed"<<std::endl;
                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif
                            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int k=(m_uiElementOrder-paddWidth);k<(m_uiElementOrder+1);k++)
                                        unzippedVec[offset+(ek*m_uiElementOrder+k-(m_uiElementOrder-paddWidth))*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];
                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c3.stop();
                            #endif

                        }

                    }

                }

                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_face[4].stop();
                #endif
                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_face[5].start();
                #endif
                
                if((pNodes[elem].maxZ()==blkNode.maxZ()))
                {
                    assert(ek==(1u<<(regLev-blkNode.getLevel()))-1);
                    lookUp=m_uiE2EMapping[elem*m_uiNumDirections+OCT_DIR_FRONT];
                    if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                    {
                        if(pNodes[lookUp].getLevel()==pNodes[elem].getLevel())
                        {

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c1.start();
                            #endif
                            assert(paddWidth<(m_uiElementOrder+1));
                            this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),lookUp);

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif
                            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int k=0;k<(paddWidth+1);k++)
                                        unzippedVec[offset+((ek+1)*m_uiElementOrder+paddWidth+k)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=lookUpElementVec[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c1.stop();
                            #endif


                        }else if(pNodes[lookUp].getLevel()<pNodes[elem].getLevel())
                        {

                            assert(pNodes[lookUp].getLevel()+1==regLev);
                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c2.start();
                            #endif
                            mid_bit=m_uiMaxDepth - pNodes[lookUp].getLevel()-1;
                            cnum=( (((((pNodes[elem].getZ()+sz)) >> mid_bit) & 1u) << 2u) | ((((pNodes[elem].getY()) >> mid_bit) & 1u) << 1u) | (((pNodes[elem].getX()) >>mid_bit) & 1u));
                            //std::cout<<"elem: "<<elem<<" : "<<m_uiAllElements[elem]<<" lookup: "<<m_uiAllElements[lookUp]<<" child: "<<ot::TreeNode(pNodes[elem].getX()+sz,pNodes[elem].getY(),pNodes[elem].getZ(),pNodes[elem].getLevel(),m_uiDim,m_uiMaxDepth)<<" cnum: "<<cnum<<std::endl;

                            

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.start();
                            #endif

                            

                            #ifdef USE_FD_INTERP_FOR_UNZIP
                            const int st = this->getBlkBdyParentCNums(blk,elem,OCT_DIR_FRONT,child.data(),fid,cid);
                            if(st > 0)
                            {
                                const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
                                this->getBlkBoundaryParentNodes(zippedVec, lookUpElementVec.data(), interpolationInput.data(), interpOrInjectionOut.data(), lookUp, fid, cid,child.data());
                                for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
                                {
                                    assert(pNodes[lookUp] == pNodes[m_uiE2EMapping[child[fid[w]]*m_uiNumDirections + OCT_DIR_FRONT]]);
                                    assert(child[fid[w]] != LOOK_UP_TABLE_DEFAULT);

                                    if(child[fid[w]]<m_uiElementLocalBegin || child[fid[w]]>=m_uiElementLocalEnd)
                                        continue;
                                    
                                    this->parent2ChildInterpolation(lookUpElementVec.data(),interpOrInjectionOut.data(),cid[w],m_uiDim);
                                    
                                    const ot::Block blk_fd = m_uiLocalBlockList[m_uiE2BlkMap[(child[fid[w]] - m_uiElementLocalBegin)]];
                                    const ot::TreeNode blkNode_fd = blk_fd.getBlockNode();
                                    const unsigned int regL_fd = blk_fd.getRegularGridLev();
                                    
                                    const unsigned int lx_fd = blk_fd.getAllocationSzX();
                                    const unsigned int ly_fd = blk_fd.getAllocationSzY();
                                    const unsigned int lz_fd = blk_fd.getAllocationSzZ();

                                    const unsigned int offset_fd = blk_fd.getOffset();

                                    const unsigned int ei_fd = (pNodes[child[fid[w]]].getX()-blkNode_fd.getX())>>(m_uiMaxDepth-regL_fd);
                                    const unsigned int ej_fd = (pNodes[child[fid[w]]].getY()-blkNode_fd.getY())>>(m_uiMaxDepth-regL_fd);
                                    const unsigned int ek_fd = (pNodes[child[fid[w]]].getZ()-blkNode_fd.getZ())>>(m_uiMaxDepth-regL_fd);

                                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int k=0;k<(paddWidth+1);k++)
                                        unzippedVec[offset_fd+((ek_fd+1)*m_uiElementOrder+paddWidth+k)*(ly_fd*lx_fd)+(ej_fd*m_uiElementOrder+j+paddWidth)*(lx_fd)+(ei_fd*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                                }

                                
                            }
                            #else
                            this->getElementNodalValues(zippedVec,&(*(lookUpElementVec.begin())),lookUp);
                            this->parent2ChildInterpolation(&(*(lookUpElementVec.begin())),&(*(interpOrInjectionOut.begin())),cnum);
                            assert(paddWidth<(m_uiElementOrder+1));

                            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int k=0;k<(paddWidth+1);k++)
                                        unzippedVec[offset+((ek+1)*m_uiElementOrder+paddWidth+k)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #endif

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c2.stop();
                            #endif



                        }else if(pNodes[lookUp].getLevel()>pNodes[elem].getLevel())
                        {
                            // get the immediate neighbours. These cannot be LOOK_UP_TABLE_DEFAULT.
                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_f_c3.start();
                            #endif
                            child[0]=lookUp;
                            child[1]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_RIGHT];
                            assert(child[1]!=LOOK_UP_TABLE_DEFAULT);
                            child[2]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_UP];
                            assert(child[2]!=LOOK_UP_TABLE_DEFAULT);
                            child[3]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_UP];
                            assert(child[3]!=LOOK_UP_TABLE_DEFAULT);

                            if(m_uiElementOrder ==4 && paddWidth==3)
                            {
                                child[4]=m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_FRONT];
                                child[5]=m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_FRONT];
                                child[6]=m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_FRONT];
                                child[7]=m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_FRONT];

                            }else
                            {
                                child[4]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[0]*m_uiNumDirections+OCT_DIR_FRONT];
                                child[5]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[1]*m_uiNumDirections+OCT_DIR_FRONT];
                                child[6]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[2]*m_uiNumDirections+OCT_DIR_FRONT];
                                child[7]=LOOK_UP_TABLE_DEFAULT;//m_uiE2EMapping[child[3]*m_uiNumDirections+OCT_DIR_FRONT];

                            }

                            this->child2ParentInjection(zippedVec,interpOrInjectionOut.data(),child.data(),pNodes[lookUp].getLevel());


                            #ifdef DEBUG_UNZIP_OP_3PT
                                faceNeighCnum1[0]=0;faceNeighCnum1[1]=1;faceNeighCnum1[2]=2;faceNeighCnum1[3]=3;
                                faceNeighCnum2[0]=4;faceNeighCnum2[1]=5;faceNeighCnum2[2]=6;faceNeighCnum2[3]=7;

                                for(unsigned int index=0;index<(NUM_CHILDREN>>1u);index++)
                                {
                                interpUpWind(fd::D1_ORDER_4_UPWIND,elem,child[faceNeighCnum1[index]],&(*(lookUpElementVec.begin())),faceNeighCnum2[index],&(*(parentEleInterpIn.begin())),&(*(parentEleInterpOut.begin())),OCT_DIR_FRONT,paddWidth,zippedVec,&(*(interpOrInjectionOut.begin())));
                                }
                            #endif

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                            dendro::timer::t_unzip_sync_cpy.start();
                            #endif
                            for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                                for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                                    for(unsigned int k=0;k<(paddWidth+1);k++)
                                        unzippedVec[offset+((ek+1)*m_uiElementOrder+paddWidth+k)*(ly*lx)+(ej*m_uiElementOrder+j+paddWidth)*(lx)+(ei*m_uiElementOrder+i+paddWidth)]=interpOrInjectionOut[k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i];

                            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                                dendro::timer::t_unzip_sync_cpy.stop();
                                dendro::timer::t_unzip_sync_f_c3.stop();
                            #endif


                        }

                    }

                }

                #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                    dendro::timer::t_unzip_sync_face[5].stop();
                #endif


            }






            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                dendro::timer::t_unzip_sync_edge.start();
            #endif
                blockDiagonalUnZip(m_uiLocalBlockList[blk],zippedVec,unzippedVec);

            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                dendro::timer::t_unzip_sync_edge.stop();
            #endif

            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                dendro::timer::t_unzip_sync_vtex.start();
            #endif
                blockVertexUnZip(m_uiLocalBlockList[blk],zippedVec,unzippedVec);

            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
                dendro::timer::t_unzip_sync_vtex.stop();
            #endif

        }

        if(m_uiElementOrder ==4 && paddWidth==3)
        {
            //std::cout<<"read spt points : "<<m_uiElementOrder<<" pwidth : "<<paddWidth<<std::endl;
            std::vector<T> recv_buf;
            recv_buf.resize(m_uiRecvOffsetRePt[m_uiActiveNpes-1] + m_uiRecvCountRePt[m_uiActiveNpes-1]);

            readSpecialPtsEnd(zippedVec,&(*(recv_buf.begin())));
            int rCount=0;
            
            for(unsigned int i=0;i<m_uiUnzip_3pt_keys.size();i++)
            {
                const std::vector<unsigned int> * ownerList = m_uiUnzip_3pt_keys[i].getOwnerList();
                for(unsigned int w=0; w< ownerList->size(); w++)
                {
                    
                    #ifdef DEBUG_UNZIP_OP_3PT
                        if(fabs(unzippedVec[(*(ownerList))[w]]-recv_buf[rCount])>1e-3)
                        {
                            std::cout<<"rank: "<<m_uiActiveRank<<" interp_deriv : "<< unzippedVec[(*(ownerList))[w]]<<" recv: "<<recv_buf[rCount]<<" diff: "<<fabs(unzippedVec[(*(ownerList))[w]]-recv_buf[rCount])<<" unzip index: "<<(*(ownerList))[w]<<std::endl;
                            //MPI_Abort(m_uiCommActive,0);
                        }
                            
                    #endif
                    unzippedVec[(*(ownerList))[w]] = recv_buf[rCount];
                }

                if(m_uiUnzip_3pt_keys[i].getOwnerList()->size()>0)
                    rCount++;

                    
            }

        }
        

        

    }

    template<typename T>
    void Mesh::readSpecialPtsBegin(const T* in)
    {

        if(m_uiGlobalNpes==1)
            return;


         // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        std::vector<T> eVec;
        eVec.resize(m_uiNpE);

        if(m_uiIsActive)
        {
            const unsigned int sendBSz=m_uiSendOffsetRePt[m_uiActiveNpes-1] + m_uiSendCountRePt[m_uiActiveNpes-1];
            const unsigned int recvBSz=m_uiRecvOffsetRePt[m_uiActiveNpes-1] + m_uiRecvCountRePt[m_uiActiveNpes-1];

            AsyncExchangeContex ctx(in);
            MPI_Comm commActive= m_uiCommActive;
            unsigned int proc_id;

            if(recvBSz)
            {
                ctx.allocateRecvBuffer((sizeof(T)*recvBSz));
                recvB=(T*)ctx.getRecvBuffer();

                // active recv procs
                for(unsigned int recv_p=0;recv_p<m_uiReqRecvProcList.size();recv_p++)
                {
                    proc_id=m_uiReqRecvProcList[recv_p];
                    MPI_Request* req=new MPI_Request();
                    par::Mpi_Irecv((recvB+m_uiRecvOffsetRePt[proc_id]),m_uiRecvCountRePt[proc_id],proc_id,m_uiCommTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }

            }

            if(sendBSz)
            {
                ctx.allocateSendBuffer(sizeof(T)*sendBSz);
                sendB=(T*)ctx.getSendBuffer();


                const unsigned int nx = m_uiElementOrder + 1;
                const unsigned int ny = m_uiElementOrder + 1;
                const unsigned int nz = m_uiElementOrder + 1;

                std::vector<unsigned int>* ownerList;
                unsigned int ownerID, ii_x, jj_y, kk_z;

                
                for(unsigned int i=0; i< m_uiUnzip_3pt_ele.size(); i++)
                {
                    ot::Key tmpEleKey= m_uiUnzip_3pt_ele[i];
                    assert((tmpEleKey.getFlag() & OCT_FOUND));
                    const unsigned int eleID = tmpEleKey.getSearchResult();
                    this->getElementNodalValues(in,&(*(eVec.begin())),eleID);
                    
                    const unsigned int step_sz = ((1u<< (m_uiMaxDepth - m_uiAllElements[eleID].getLevel()))/m_uiElementOrder);
                    ownerList = tmpEleKey.getOwnerList();
                    for(unsigned int w=0; w< ownerList->size();w++)
                    {
                        
                        const unsigned int ii = (m_uiUnzip_3pt_recv_keys[(*ownerList)[w]].minX() - m_uiAllElements[eleID].minX())/(step_sz); 
                        const unsigned int jj = (m_uiUnzip_3pt_recv_keys[(*ownerList)[w]].minY() - m_uiAllElements[eleID].minY())/(step_sz); 
                        const unsigned int kk = (m_uiUnzip_3pt_recv_keys[(*ownerList)[w]].minZ() - m_uiAllElements[eleID].minZ())/(step_sz);
                        const std::vector<unsigned int > * ownerList1 = m_uiUnzip_3pt_recv_keys[(*ownerList)[w]].getOwnerList();

                        for(unsigned int w1 = 0; w1 < ownerList1->size() ; w1++)
                        {
                            // if(m_uiActiveRank==1 && (*ownerList1)[w1]<18 )
                            //     std::cout<<" rank: "<<m_uiActiveRank<<" putting : "<<m_uiUnzip_3pt_recv_keys[(*ownerList)[w]]<<" to send buf loc: "<<(*ownerList1)[w1]<<std::endl;

                            sendB[(*ownerList1)[w1]] = eVec[kk * ny * nx + jj * nx + ii];
                        }
                            
                        
                    }
                }


                // active send procs
                for(unsigned int send_p=0;send_p<m_uiReqSendProcList.size();send_p++)
                {
                    proc_id=m_uiReqSendProcList[send_p];
                    MPI_Request * req=new MPI_Request();
                    par::Mpi_Isend(sendB+m_uiSendOffsetRePt[proc_id],m_uiSendCountRePt[proc_id],proc_id,m_uiCommTag,commActive,req);
                    ctx.getRequestList().push_back(req);

                }


            }

            m_uiCommTag++;
            m_uiMPIContexts.push_back(ctx);



        }

    }

    template <typename T>
    void Mesh::readSpecialPtsEnd(const T *in, T* out)
    {
        if(m_uiGlobalNpes == 1)
            return;

        // send recv buffers.
        T* sendB = NULL;
        T* recvB = NULL;

        if(m_uiIsActive)
        {
            const unsigned int sendBSz=m_uiSendOffsetRePt[m_uiActiveNpes-1] + m_uiSendCountRePt[m_uiActiveNpes-1];
            const unsigned int recvBSz=m_uiRecvOffsetRePt[m_uiActiveNpes-1] + m_uiRecvCountRePt[m_uiActiveNpes-1];

            //std::cout<<"rank: "<<m_uiActiveRank<<" recv sz: "<<recvBSz<<std::endl;

            unsigned int proc_id;
            unsigned int ctxIndex=0;

            for(unsigned int i=0;i<m_uiMPIContexts.size();i++)
            {
                if(m_uiMPIContexts[i].getBuffer()==in)
                {
                    ctxIndex=i;
                    break;
                }

            }

            MPI_Status status;
            // need to wait for the commns to finish ...
            for (unsigned int i = 0; i < m_uiMPIContexts[ctxIndex].getRequestList().size(); i++) {
                MPI_Wait(m_uiMPIContexts[ctxIndex].getRequestList()[i], &status);
            }

            if(recvBSz)
            {
                // copy the recv data to the vec
                recvB=(T*)m_uiMPIContexts[ctxIndex].getRecvBuffer();
                std::memcpy(out,recvB,sizeof(T)*recvBSz);
                
                // for(unsigned int i=0; i<recvBSz ; i++ )
                //     out[i] = recvB[i];
            }



            m_uiMPIContexts[ctxIndex].deAllocateSendBuffer();
            m_uiMPIContexts[ctxIndex].deAllocateRecvBuffer();

            for (unsigned int i = 0; i < m_uiMPIContexts[ctxIndex].getRequestList().size(); i++)
                delete m_uiMPIContexts[ctxIndex].getRequestList()[i];

            m_uiMPIContexts[ctxIndex].getRequestList().clear();

            // remove the context ...
            m_uiMPIContexts.erase(m_uiMPIContexts.begin() + ctxIndex);


        }

        return;
    }


    

    template<typename T>
    int Mesh::getFaceNeighborValues(unsigned int eleID, const T* in, T* out, T* coords, unsigned int * neighID, unsigned int face, NeighbourLevel & level) const
    {

        if(!m_uiIsActive)
            return (0);

        const unsigned int lookUp=m_uiE2EMapping[eleID*m_uiNumDirections + face];
        if(lookUp==LOOK_UP_TABLE_DEFAULT)
            return (0);

        const unsigned int l1=m_uiAllElements[eleID].getLevel();
        const unsigned int l2=m_uiAllElements[lookUp].getLevel();

        for(unsigned int i=0;i<(NUM_CHILDREN>>1);i++)
            neighID[i]=LOOK_UP_TABLE_DEFAULT;

        int num_face_neighbours = 1;
        if(l1==l2)
        {
            // both elements are in the same level.
            level = NeighbourLevel::SAME;
            neighID[0]=lookUp;
            this->getElementNodalValues(in,out,lookUp);

            // coordinate computation
            const ot::TreeNode lookUpOct=m_uiAllElements[lookUp];
            const unsigned int sz=1u<<(m_uiMaxDepth-lookUpOct.getLevel());

            for(unsigned int k=0;k< (m_uiElementOrder+1);k++ )
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    {
                        coords[m_uiDim* (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 0] = lookUpOct.minX() + i * (sz/(T)m_uiElementOrder);
                        coords[m_uiDim* (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 1] = lookUpOct.minY() + j * (sz/(T)m_uiElementOrder);
                        coords[m_uiDim* (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 2] = lookUpOct.minZ() + k * (sz/(T)m_uiElementOrder);
                    }
        }else if(l2<l1)
        {
            level = NeighbourLevel::COARSE;
            // lookUp octant is coaser than eleID.
            neighID[0]=lookUp;
            this->getElementNodalValues(in,out + m_uiNpE ,lookUp);

            unsigned int sz1=1u<<(m_uiMaxDepth-m_uiAllElements[eleID].getLevel());

            unsigned int x = m_uiAllElements[eleID].minX();
            unsigned int y = m_uiAllElements[eleID].minY();
            unsigned int z = m_uiAllElements[eleID].minZ();



            ot::TreeNode tmpOct;
            unsigned int cnum;
            switch (face)
            {
                case OCT_DIR_LEFT:
                    tmpOct=ot::TreeNode(x-sz1,y,z,l1,m_uiDim,m_uiMaxDepth);
                    cnum=tmpOct.getMortonIndex();
                    break;

                case OCT_DIR_RIGHT:
                    tmpOct=ot::TreeNode(x+sz1,y,z,l1,m_uiDim,m_uiMaxDepth);
                    cnum=tmpOct.getMortonIndex();
                    break;

                case OCT_DIR_DOWN:
                    tmpOct=ot::TreeNode(x,y-sz1,z,l1,m_uiDim,m_uiMaxDepth);
                    cnum=tmpOct.getMortonIndex();
                    break;

                case OCT_DIR_UP:
                    tmpOct=ot::TreeNode(x,y+sz1,z,l1,m_uiDim,m_uiMaxDepth);
                    cnum=tmpOct.getMortonIndex();
                    break;


                case OCT_DIR_BACK:
                    tmpOct=ot::TreeNode(x,y,z-sz1,l1,m_uiDim,m_uiMaxDepth);
                    cnum=tmpOct.getMortonIndex();
                    break;

                case OCT_DIR_FRONT:
                    tmpOct=ot::TreeNode(x,y,z+sz1,l1,m_uiDim,m_uiMaxDepth);
                    cnum=tmpOct.getMortonIndex();
                    break;

                default:
                    std::cout<<"global rank : "<<m_uiGlobalRank<<" dir: "<<face<<" is invalid. Function : "<<__func__<<std::endl;
                    MPI_Abort(m_uiCommGlobal,0);
                    break;
            }


            this->parent2ChildInterpolation(out+m_uiNpE,out,cnum,m_uiDim);

            // coordinate computation
            const ot::TreeNode lookUpOct=tmpOct;
            const unsigned int sz=1u<<(m_uiMaxDepth-lookUpOct.getLevel());

            for(unsigned int k=0;k< (m_uiElementOrder+1);k++ )
                for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                    for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                    {
                        coords[m_uiDim* (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 0] = lookUpOct.minX() + i * (sz/(T)m_uiElementOrder);
                        coords[m_uiDim* (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 1] = lookUpOct.minY() + j * (sz/(T)m_uiElementOrder);
                        coords[m_uiDim* (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 2] = lookUpOct.minZ() + k * (sz/(T)m_uiElementOrder);
                    }


        }else{
            // lookUp octant is finer than eleID.
            assert(l2>l1);

            unsigned int dir, dirOp, dir1,dir2;
            num_face_neighbours = 4;
            level = NeighbourLevel::REFINE;
            switch (face)
            {

                case OCT_DIR_LEFT:

                    dir=OCT_DIR_LEFT;
                    dirOp=OCT_DIR_RIGHT;
                    dir1=OCT_DIR_FRONT;
                    dir2=OCT_DIR_UP;

                    neighID[0]=lookUp;
                    neighID[1]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir2];
                    neighID[2]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir1];
                    neighID[3]=m_uiE2EMapping[neighID[1]*NUM_FACES+dir1];


                    break;

                case OCT_DIR_RIGHT:

                    dir=OCT_DIR_RIGHT;
                    dirOp=OCT_DIR_LEFT;

                    dir1=OCT_DIR_FRONT;
                    dir2=OCT_DIR_UP;

                    neighID[0]=lookUp;
                    neighID[1]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir2];
                    neighID[2]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir1];
                    neighID[3]=m_uiE2EMapping[neighID[1]*NUM_FACES+dir1];

                    break;

                case OCT_DIR_DOWN:

                    dir=OCT_DIR_DOWN;
                    dirOp=OCT_DIR_UP;

                    dir1=OCT_DIR_FRONT;
                    dir2=OCT_DIR_RIGHT;


                    neighID[0]=lookUp;
                    neighID[1]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir2];
                    neighID[2]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir1];
                    neighID[3]=m_uiE2EMapping[neighID[1]*NUM_FACES+dir1];

                    break;

                case OCT_DIR_UP:

                    dir=OCT_DIR_UP;
                    dirOp=OCT_DIR_DOWN;

                    dir1=OCT_DIR_FRONT;
                    dir2=OCT_DIR_RIGHT;

                    neighID[0]=lookUp;
                    neighID[1]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir2];
                    neighID[2]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir1];
                    neighID[3]=m_uiE2EMapping[neighID[1]*NUM_FACES+dir1];


                    break;

                case OCT_DIR_BACK:

                    dir=OCT_DIR_BACK;
                    dirOp=OCT_DIR_FRONT;

                    dir1=OCT_DIR_UP;
                    dir2=OCT_DIR_RIGHT;

                    neighID[0]=lookUp;
                    neighID[1]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir2];
                    neighID[2]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir1];
                    neighID[3]=m_uiE2EMapping[neighID[1]*NUM_FACES+dir1];


                    break;

                case OCT_DIR_FRONT:

                    dir=OCT_DIR_FRONT;
                    dirOp=OCT_DIR_BACK;

                    dir1=OCT_DIR_UP;
                    dir2=OCT_DIR_RIGHT;


                    neighID[0]=lookUp;
                    neighID[1]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir2];
                    neighID[2]=m_uiE2EMapping[neighID[0]*NUM_FACES+dir1];
                    neighID[3]=m_uiE2EMapping[neighID[1]*NUM_FACES+dir1];


                    break;

                default:
                    std::cout<<"global rank : "<<m_uiGlobalRank<<" dir: "<<face<<" is invalid. Function : "<<__func__<<std::endl;
                    MPI_Abort(m_uiCommGlobal,0);


            }

            for(unsigned int child=0;child<(NUM_CHILDREN>>1);child++)
            {

                this->getElementNodalValues(in,out+child*m_uiNpE,neighID[child]);

                // coordinate computation
                const ot::TreeNode lookUpOct=m_uiAllElements[neighID[child]];
                const unsigned int sz=1u<<(m_uiMaxDepth-lookUpOct.getLevel());

                for(unsigned int k=0;k< (m_uiElementOrder+1);k++ )
                    for(unsigned int j=0;j<(m_uiElementOrder+1);j++)
                        for(unsigned int i=0;i<(m_uiElementOrder+1);i++)
                        {
                            coords[ child*m_uiNpE*m_uiDim + m_uiDim * (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 0] = lookUpOct.minX() + i * (sz/(T)m_uiElementOrder);
                            coords[ child*m_uiNpE*m_uiDim + m_uiDim * (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 1] = lookUpOct.minY() + j * (sz/(T)m_uiElementOrder);
                            coords[ child*m_uiNpE*m_uiDim + m_uiDim * (k*(m_uiElementOrder+1)*(m_uiElementOrder+1)+j*(m_uiElementOrder+1)+i) + 2] = lookUpOct.minZ() + k * (sz/(T)m_uiElementOrder);
                        }
            }




        }

        return num_face_neighbours;

    }

    template<typename T>
    void Mesh::getUnzipElementalNodalValues(const T* uzipVec, unsigned int blkID, unsigned int ele, T*out, bool isPadded) const
    {
        const ot::Block block = m_uiLocalBlockList[blkID];
        ot::TreeNode blkNode = m_uiLocalBlockList[blkID].getBlockNode();
        const unsigned int eleBegin = block.getLocalElementBegin();
        const unsigned int eleEnd = block.getLocalElementEnd();

        assert(eleBegin <= ele && ele < eleEnd);
        const unsigned int regLev=block.getRegularGridLev();
        const unsigned int lx=block.getAllocationSzX();
        const unsigned int ly=block.getAllocationSzY();
        const unsigned int lz=block.getAllocationSzZ();
        const unsigned int offset=block.getOffset();
        const unsigned int paddWidth=block.get1DPadWidth();

        const unsigned int ei=(m_uiAllElements[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
        const unsigned int ej=(m_uiAllElements[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
        const unsigned int ek=(m_uiAllElements[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

        

        if(isPadded)
        {
            const unsigned int ib = ei*m_uiElementOrder;
            const unsigned int ie = ei*m_uiElementOrder + (m_uiElementOrder+1) + 2*paddWidth ;

            const unsigned int jb = ej*m_uiElementOrder;
            const unsigned int je = ej*m_uiElementOrder + (m_uiElementOrder+1) + 2*paddWidth ;

            const unsigned int kb = ek*m_uiElementOrder;
            const unsigned int ke = ek*m_uiElementOrder + (m_uiElementOrder+1) + 2*paddWidth ;

            const unsigned int en[3] = {(m_uiElementOrder+1) + 2*paddWidth , (m_uiElementOrder+1) + 2*paddWidth, (m_uiElementOrder+1) + 2*paddWidth }; 

            for(unsigned int k=kb; k< ke; k++)
             for(unsigned int j=jb; j< je; j++)
              for(unsigned int i=ib; i< ie; i++)
               out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[offset + k*ly*lx + j*lx + i];

        }else
        {
            const unsigned int ib = ei*m_uiElementOrder + paddWidth;
            const unsigned int ie = ei*m_uiElementOrder + (m_uiElementOrder+1) ;

            const unsigned int jb = ej*m_uiElementOrder + paddWidth ;
            const unsigned int je = ej*m_uiElementOrder + (m_uiElementOrder+1) ;

            const unsigned int kb = ek*m_uiElementOrder + paddWidth ;
            const unsigned int ke = ek*m_uiElementOrder + (m_uiElementOrder+1) ;

            const unsigned int en[3] = {(m_uiElementOrder+1) , (m_uiElementOrder+1) , (m_uiElementOrder+1) };

            for(unsigned int k=kb; k< ke; k++)
             for(unsigned int j=jb; j< je; j++)
              for(unsigned int i=ib; i< ie; i++)
                out[(k-kb) * en[1]*en[0] + (j-jb)*en[1] + (i-ib)] = uzipVec[offset + k*ly*lx + j*lx + i];



        }
        
        
    }


    template<typename T>
    void Mesh::getBlkBoundaryParentNodes(const T* zipVec, T* out, T*w1, T*w2, unsigned int lookUp, const unsigned int * fid, const unsigned int* cid,const unsigned int * child)
    {

        const unsigned int NUM_CHILDREN_BY2 = (NUM_CHILDREN>>1u);
        const unsigned int eorder_by2 = (m_uiElementOrder+1)>>1u;
        const unsigned int nx = m_uiElementOrder + 1;
        const unsigned int ny = m_uiElementOrder + 1;
        const unsigned int nz = m_uiElementOrder + 1;
        
        unsigned char bit[3];

        // finner elements. 
        for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
        {
            this->getElementNodalValues(zipVec,w1,child[fid[w]]);
            //std::cout<<" cnum : "<<fid[w]<<std::endl;
            bit[0] = binOp::getBit(fid[w],0);
            bit[1] = binOp::getBit(fid[w],1);
            bit[2] = binOp::getBit(fid[w],2);

            const unsigned int kb = bit[2] * eorder_by2 ;
            unsigned int ke = kb + eorder_by2 +1;

            const unsigned int jb = bit[1] * eorder_by2 ;
            unsigned int je = jb + eorder_by2 +1;

            const unsigned int ib = bit[0] * eorder_by2 ;
            unsigned int ie = ib + eorder_by2 +1;

            for(unsigned int k=0; k< nz; k+=2)
                for(unsigned int j=0; j < ny; j+=2)
                for(unsigned int i=0; i < nx; i+=2)
                out[(kb + (k>>1u))*ny*nx + (jb + (j>>1u))*nx + (ib + (i>>1u))] = w1[k*ny*nx + j*nx + i];

        }

        this->getElementNodalValues(zipVec,w1,lookUp);
        // coarser elements. 
        for(unsigned int w =0; w < NUM_CHILDREN_BY2 ; w++)
        {
            this->parent2ChildInterpolation(w1,w2,fid[w],m_uiDim);

            //std::cout<<" cnum : "<<fid[w]<<std::endl;
            bit[0] = binOp::getBit(cid[w],0);
            bit[1] = binOp::getBit(cid[w],1);
            bit[2] = binOp::getBit(cid[w],2);

            const unsigned int kb = bit[2] * eorder_by2 ;
            unsigned int ke = kb + eorder_by2 +1;

            const unsigned int jb = bit[1] * eorder_by2 ;
            unsigned int je = jb + eorder_by2 +1;

            const unsigned int ib = bit[0] * eorder_by2 ;
            unsigned int ie = ib + eorder_by2 +1;

            for(unsigned int k=0; k< nz; k+=2)
                for(unsigned int j=0; j < ny; j+=2)
                for(unsigned int i=0; i < nx; i+=2)
                {
                    // std::cout<< " cnum : "<<fid[w]<< " left lookup value: ijk: "<<i<<j<<k<<" "<<lookUpElementVec[(kb + (k>>1u))*ny*nx + (jb + (j>>1u))*nx + (ib + (i>>1u))]<< " inject value: "<<w2[k*ny*nx + j*nx + i]<<"";
                    // printf("lookup idx (%d,%d,%d)\n",(ib + (i>>1u)),(jb + (j>>1u)), (kb + (k>>1u)));
                    out[(kb + (k>>1u))*ny*nx + (jb + (j>>1u))*nx + (ib + (i>>1u))] = w2[k*ny*nx + j*nx + i];
                }
                

        }

    }

}//end of namespase ot





