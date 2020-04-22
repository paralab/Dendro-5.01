/**
 * @file blkAsync.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Fully asynchronous storage to ot::Block data structure for non-uniform time stepping. 
 * @version 0.1
 * @date 2020-03-17
 * 
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2020
 * 
 */

#pragma once
#include "ts.h"
#include "mesh.h"
#include "enutsOp.h"

namespace ts
{
    enum BLK_ASYNC_VEC_MODE{BLK_CG=0, BLK_DG, BLK_UNZIP};


    template<typename T>
    class BlockAsyncVector
    {

        protected: 

            /**@brief: block ID (local)*/
            unsigned int m_uiBlkID = LOOK_UP_TABLE_DEFAULT;

            /**@brief : data vector (rhs in)*/
            T* m_uiVec = NULL;

            unsigned int m_uiSz[3];

            /**@brief: is block synced and ready to be evolved. */
            bool m_uiIsSynced = false;

            /**@brief: type of the vector. */
            BLK_ASYNC_VEC_MODE m_uiMode;

            /**@brief: number of dof. */
            unsigned int m_uiDof;

            /**@brief: list of send requests*/
            std::vector<MPI_Request*> m_uiSendReq;

            /**@brief: list of recv requests*/
            //std::vector<MPI_Request*> m_uiRecvReq;

            
            
        public:
            
            /**@brief: default constructor.*/
            BlockAsyncVector(){};

            /**@brief: default destructor*/
            ~BlockAsyncVector()
            {
                if(m_uiVec!=NULL);
                    delete [] m_uiVec;

                m_uiVec==NULL;
                
            };
            
            /**
             * @brief allocate memory for a block vector
             * @param blk :block id
             * @param sz  : size of the block x,y,z
             * @param synced : true if the block vector is synced 
             * @param mode : BLK_ASYNC_VEC_MODE 
             * @param dof : number of dof
             */
            void createVec(unsigned int blk, const unsigned int* sz, bool synced=false, BLK_ASYNC_VEC_MODE mode = BLK_ASYNC_VEC_MODE::BLK_UNZIP, unsigned int dof=1)
            {
                m_uiBlkID = blk ;
                m_uiSz[0] = sz[0]; m_uiSz[1] = sz[1]; m_uiSz[2] = sz[2]; 
                const unsigned int NN = (m_uiSz[0]*m_uiSz[1]*m_uiSz[2]);
                m_uiIsSynced = synced;
                m_uiMode = mode;
                m_uiDof = dof;

                if(NN==0)
                {
                    m_uiVec=NULL;
                    return;
                }

                m_uiVec = new T [m_uiDof * NN ];

                return;

            }

            /**
             * @brief : destroy a block vector. 
             */
            void destroyVec()
            {
                m_uiBlkID = LOOK_UP_TABLE_DEFAULT;
                m_uiIsSynced = false;
        
                delete [] m_uiVec;
                m_uiVec=NULL;

            }

            /**
             * @brief copy the block buffer from an unzip vector. 
             * @param pMesh : mesh data strucuture. 
             * @param vUnzip : unzip vector. 
             * @param isWithPadding : if true the copy will happen with padding. otherwise padding region is not coppied. 
             * @param dof : number of dofs
             */
            void copyFromUnzip(const ot::Mesh*pMesh, const T* const vUnzip , bool isWithPadding, unsigned int dof=1)
            {
                if(!(pMesh->isActive()))
                    return;
                
                const unsigned int blk = m_uiBlkID;
                const ot::Block* blkList  = pMesh->getLocalBlockList().data();
                
                if(m_uiMode == BLK_ASYNC_VEC_MODE::BLK_UNZIP)
                {
                    const unsigned int unzipSz = pMesh->getDegOfFreedomUnZip();
                    const unsigned int offset  = blkList[blk].getOffset();
                    const unsigned int nx      = blkList[blk].getAllocationSzX();
                    const unsigned int ny      = blkList[blk].getAllocationSzY();
                    const unsigned int nz      = blkList[blk].getAllocationSzZ();
                    const unsigned int pw      = blkList[blk].get1DPadWidth();

                    assert(m_uiSz[0]==nx && m_uiSz[1] ==ny && m_uiSz[2]==nz);

                    if(isWithPadding)
                    {
                        for(unsigned int v =0; v < dof; v++ )
                        {
                            const T* v1 = vUnzip + v * unzipSz;
                            std::memcpy(m_uiVec + v*nx*ny*nz ,&v1[offset],sizeof(T)*nx*ny*nz);
                        }


                    }else
                    {

                        for(unsigned int v =0; v < dof; v++ )
                        {
                            const T* v1 = vUnzip + v * unzipSz;
                            for(unsigned int k= pw; k < nz-pw; k++)
                            for(unsigned int j= pw; j < ny-pw; j++)
                            for(unsigned int i= pw; i < nx-pw; i++)
                                m_uiVec[(v*nx*ny*nz) + k*ny*nx + j*nx + i] = v1[offset + k * ny*nx +  j * nx + i];
                        }
                        
                    }
                    


                }

                return ;
            }

            /**
             * @brief copy a given other vector to the this vector. 
             * @param other : BlockAsyncVec that need to be coppied. 
             */
            void vecCopy(const BlockAsyncVector<T>& other)
            {
                m_uiBlkID = other.m_uiBlkID;
                m_uiSendReq = other.m_uiSendReq;
                m_uiIsSynced = other.isSynced();
                m_uiMode = other.m_uiMode;
                m_uiDof = other.m_uiDof;

                m_uiSz[0] = other.m_uiSz[0];
                m_uiSz[1] = other.m_uiSz[1];
                m_uiSz[2] = other.m_uiSz[2];

                const unsigned int NN = m_uiSz[0]*m_uiSz[1]*m_uiSz[2];
                
                if(m_uiMode == BLK_ASYNC_VEC_MODE::BLK_UNZIP)
                    std::memcpy(m_uiVec,other.m_uiVec,sizeof(T)*NN*m_uiDof);

                return;
            }


            /**
             * @brief perform zip operation for the block
             * @param pMesh : octree Mesh object.
             * @param zipVec : zipped vector
             * @param dof : number of dof. 
             */
            void zip(ot::Mesh*pMesh, T* zipVec,  unsigned int dof=1) const
            {

                if(!(pMesh->isActive()))
                    return;

                const ot::Block* blkList = pMesh->getLocalBlockList().data();
                const unsigned int * e2n = pMesh->getE2NMapping().data();
                const unsigned int lx = m_uiSz[0];
                const unsigned int ly = m_uiSz[1];
                const unsigned int lz = m_uiSz[2];
                
                const unsigned int paddWidth=blkList[m_uiBlkID].get1DPadWidth();

                const unsigned int eOrder = pMesh->getElementOrder();
                const unsigned int nPe = pMesh->getNumNodesPerElement();

                const unsigned int vsz_cg = pMesh->getDegOfFreedom();
                const ot::TreeNode* pNodes = pMesh->getAllElements().data();
                

                
                const ot::TreeNode blkNode =  blkList[m_uiBlkID].getBlockNode();
                const unsigned int regLev  =  blkList[m_uiBlkID].getRegularGridLev();
                
                if(m_uiMode == BLK_ASYNC_VEC_MODE::BLK_UNZIP)
                {
                    for(unsigned int v=0; v < dof; v++)
                    {
                        
                        for(unsigned int elem = blkList[m_uiBlkID].getLocalElementBegin(); elem < blkList[m_uiBlkID].getLocalElementEnd(); elem++)
                        {
                            const unsigned int ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                            const unsigned int ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                            const unsigned int ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                            for(unsigned int k=0; k < (eOrder+1); k++)
                            for(unsigned int j=0; j < (eOrder+1); j++)
                            for(unsigned int i=0; i < (eOrder+1); i++)
                            {

                                const bool isHanging = pMesh->isNodeHanging(elem,i,j,k);
                                if(!isHanging)
                                    zipVec[ (v*vsz_cg) + e2n[elem*nPe + k*(eOrder+1)*(eOrder+1) + j* (eOrder+1)+ i]] = m_uiVec[ (v*lx*ly*lz) + (ek*eOrder+k+paddWidth)*(ly*lx)+(ej*eOrder+j+paddWidth)*(lx)+(ei*eOrder+i+paddWidth)];
                                else
                                {
                                    const unsigned int cnum=pNodes[(elem)].getMortonIndex();
                                    const unsigned int iix = eOrder * (int) (cnum & 1u)  +  i;
                                    const unsigned int jjy = eOrder * (int) ((cnum & 2u)>>1u)  +  j;
                                    const unsigned int kkz = eOrder * (int) ((cnum & 4u)>>2u)  +  k;

                                    if( (iix %2 ==0) && (jjy%2 ==0) && (kkz%2==0))
                                        zipVec[(v*vsz_cg) + e2n[elem*nPe + (kkz>>1u)*(eOrder+1)*(eOrder+1) + (jjy>>1u) * (eOrder+1)+ (iix>>1u)]] = m_uiVec[ (v*lx*ly*lz)  +   (ek*eOrder+k+paddWidth)*(ly*lx)+(ej*eOrder+j+paddWidth)*(lx)+(ei*eOrder+i+paddWidth)];

                                }
                                
                            }
            
                        }
                    
                    }
                
                }

            }

            /**
             * @brief copy block vector to a DG vector. 
             * @param pMesh 
             * @param dgVec 
             * @param dof 
             */
            void zipDG(ot::Mesh*pMesh, T* dgVec,  unsigned int dof=1) const
            {

                if(!(pMesh->isActive()))
                    return;

                const ot::Block* blkList     =   pMesh->getLocalBlockList().data();
                const unsigned int regLev    =   blkList[m_uiBlkID].getRegularGridLev();
                const ot::TreeNode* pNodes   =   pMesh->getAllElements().data();
                const ot::TreeNode blkNode   =   blkList[m_uiBlkID].getBlockNode();
                const unsigned int eOrder    =   pMesh->getElementOrder();
                const unsigned int paddWidth =   blkList[m_uiBlkID].get1DPadWidth();

                const unsigned int nPe       =  (eOrder+1)*(eOrder+1)*(eOrder+1);

                const unsigned int lx   =   blkList[m_uiBlkID].getAllocationSzX();
                const unsigned int ly   =   blkList[m_uiBlkID].getAllocationSzY();
                const unsigned int lz   =   blkList[m_uiBlkID].getAllocationSzZ();

                const unsigned int dgSz = pMesh->getAllElements().size() * pMesh->getNumNodesPerElement();
        

                for(unsigned int v=0; v < dof; v++)
                {
                    for(unsigned int elem = blkList[m_uiBlkID].getLocalElementBegin(); elem < blkList[m_uiBlkID].getLocalElementEnd(); elem++)
                    {
                        const unsigned int ei  =  (pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                        const unsigned int ej  =  (pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                        const unsigned int ek  =  (pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                        for(unsigned int k=0; k < (eOrder+1); k++)
                         for(unsigned int j=0; j < (eOrder+1); j++)
                          for(unsigned int i=0; i < (eOrder+1); i++)
                            dgVec[ (v*dgSz + elem*nPe)  +  k*(eOrder+1)*(eOrder+1) +  j*(eOrder+1) + i] = m_uiVec[ (v*lx*ly*lz)  + (ek*eOrder+k+paddWidth)*(ly*lx)+(ej*eOrder+j+paddWidth)*(lx)+(ei*eOrder+i+paddWidth)];
                    }
                
                }

                return;

            }

           

            /**@brief: returns the m_uiIsSynced status*/
            inline bool isSynced() const { return m_uiIsSynced; }
            
            /**
             * @brief compute the rhs vector based on the given application ctx. 
             * @param appCtx : const pointer to the application contex. 
             */
            template<typename ACtx>
            void computeVec(const ACtx * appCtx, const BlockAsyncVector<T>& in, T time)
            {
                const ot::Mesh* pMesh = appCtx->get_mesh();
                if(!(pMesh->isActive()))
                    return;
                
                m_uiMode = in.m_uiMode;
                m_uiDof = in.m_uiDof;
                m_uiIsSynced = false;
                ((ACtx*)appCtx)->rhs_blk(in.m_uiVec,m_uiVec,m_uiDof,m_uiBlkID,time);
                return;
            }

            /**@brief: Initiate an asynchronous send from my blocks to other processors*/
            virtual void iSend(const ot::Mesh* pMesh, T* vBuff, const  std::vector<unsigned int>& sp,  int tag, MPI_Comm comm) const 
            {

            }

            /**@brief: Wait until the send req. are completed. */
            virtual void waitOnSend()
            {

            }

            /**@brief: returns the m_uiVec pointer. */
            inline const T* data() { return m_uiVec ; }

            /**@brief: return the block id*/
            inline unsigned int getBlkID() const {return m_uiBlkID;}
            
            /**@brief: return the DOF*/
            inline unsigned int getDOF() const {return m_uiDof;}

            /**@brief: mark the block vector synced. */
            inline void mark_synced() { m_uiIsSynced =true;}

            /**@brief: mark the block vector unsynced. */
            inline void mark_unsynced() { m_uiIsSynced =false;}

            /**@brief: returns the size of the block vector. */
            inline unsigned int getSz() const { return (m_uiSz[0]*m_uiSz[1]*m_uiSz[2]); }

            /**@brief: dump out the vector. */
            void dump_vec(std::ostream & sout)
            {
                const unsigned int nn = m_uiSz[0]*m_uiSz[1]*m_uiSz[2];

                for(unsigned int v =0; v < m_uiDof; v++)
                {
                    sout<<"var : "<<v<<std::endl;
                    for(unsigned int n=0; n < nn; n++)
                        sout<<"m_uiVec["<<n<<"] : "<<m_uiVec[ v*nn + n ]<<std::endl;
                }
            }


    };





    



}// end of namespace ts.


