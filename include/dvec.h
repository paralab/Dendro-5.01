/**
 * @file dvec.h
 * @author Milinda Fernando
 * @brief Vector class with additional auxiliary information. 
 * @version 0.1
 * @date 2019-12-19
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#pragma once
#include <iostream>
#include "mpi.h"
#include "mesh.h"
#include "mathUtils.h"
namespace ot
{
    template<typename T,typename I>
    class DVector
    {   
        protected:

            /**@brief ptr for the data*/
            T* m_uiData = NULL;

            /**@brief size of the vector*/
            I m_uiSize = 0 ;

            /** @brief : true if allocated as a unzip vector format false otherwise*/
            bool m_uiIsUnzip = false;

            /**@brief : true if allocated as a elemental vector, if as an nodal vector false. */
            bool m_uiIsElemental = false;

            /**@brief : true if allocated with ghost/halo regions, false otherwise */
            bool m_uiIsGhosted = true;

            /**@brief: true if the current vector is CG vector*/
            bool m_uiIsCG =true;

            /**@brief: number of degrees of freedoms. */
            unsigned int m_uiDof;

            /**@brief: MPI Communicator for the vector*/
            MPI_Comm m_uiComm;
            
        public:
            
            /**@brief: default constructor*/
            DVector();
            
            /**@brief: destructor*/
            ~DVector();

            /**
             * @brief Construct a new vec create object
             * 
             * @param pMesh : pointer to underlying ot::Mesh object
             * @param isGhosted : true if need to create a vector with halo regions. 
             * @param isUnzip : true if this is unzip vector format (block local vector format)
             * @param isElemental : true if this is elemental vector
             * @param dof : number of degrees of freedoms. 
             */
            void VecCreate(const ot::Mesh* pMesh, bool isGhosted = false, bool isUnzip =false, bool isElemental = false, unsigned int dof=1);

            /**
             * @brief creates a element DG vector
             * 
             * @param pMesh Mesh object. 
             * @param isGhosted true if the vector needs to be ghosted DG
             * @param dof : degrees of freedom. 
             */
            void VecCreateDG(const ot::Mesh* pMesh, bool isGhosted = false , unsigned int dof=1);

            /**
             * @brief deallocates the vector object
             */
            void VecDestroy();  

            /**@brief: returns the vec pointer*/
            inline T*& GetVecArray() { return m_uiData; };
            
            /**@brief: Get the vector pointer.*/
            inline void GetVecArray(T*& vec) { vec= m_uiData; }

            /**
             * @brief update the vector pointer object. 
             * @param vec : vec pointer
             */
            inline void VecRestoreArray(T* vec) { m_uiData = vec;};

            /**
             * @brief copy vec v to current vector
             * 
             * @param v : other vector v. 
             * @param isAlloc : true is the current vector (*this) is allocated. 
             */
            void VecCopy(DVector<T,I> v, bool isAlloc = false);

            /**
             * @brief performs vector fuse multiplication and addition. 
             * 
             * @param pMesh underlying mesh data strucuture. 
             * @param v vector v
             * @param a scale parameter for aU + bV
             * @param b scale parameter for aU + bV
             * @param localOnly if true skip the ghosted regions of the array if applicable. 
             */
            void VecFMA( const ot::Mesh * pMesh , DVector<T,I> v, T a, T b, bool localOnly = true);

            /**
             * @brief 
             * @param pMesh underlying mesh data strucuture. 
             * @param a scale parameter for aU + bV
             * @param b scale parameter for aU + bV
             * @param localOnly if true skip the ghosted regions of the array if applicable.
             */
            void VecFMA( const ot::Mesh * pMesh , T a, T b, bool localOnly = true);

            /**
             * @brief : equal operator for the DVector
             * 
             * @param other : DVector object 
             * @return true if points to the same pointer. 
             * @return false otherwise. 
             */
            bool operator==(DVector<T,I>& other) const { return ( m_uiData == (other.GetVecArray())); }

            void operator=(DVector<T,I> other); 

            /**@brief: get is ghosted variable. */
            inline bool IsGhosted() const { return m_uiIsGhosted; } 

            /**@brief: true if this is an unzip vector*/
            inline bool IsUnzip() const { return m_uiIsUnzip;}

            /**@brief: true if this is an elemental vector. */
            inline bool IsElemental() const {return m_uiIsElemental;}

            /**@brief: returns the number of degrees of freedoms. */
            inline unsigned int GetDof() const {return m_uiDof;}

            inline MPI_Comm GetMPIComm() const {return m_uiComm;}

            /**
             * @brief Get the Size vector
             * @return I 
             */
            inline I GetSize() const { return m_uiSize;}

            /**@brief:
             * get the 2D array for a multiple dof vector (note that the v2d should be a pointer in heap. )
             * @param isAlloc: true if the mem allocated for T*. False otherwise. 
            */
            void Get2DArray(T** & v2d, bool isAlloc = false);

            /**
             * @brief returns the 2D vector. Note that the v2d assumeed to allocated no need to be on the heap;
             * @param v2d 
             */
            void Get2DVec(T** v2d);

            /**@brief: get size per dof*/
            inline I GetSizePerDof() const {return (m_uiSize/m_uiDof) ;}

            /**@brief: computes the min and max of the vector. */
            void VecMinMax(const ot::Mesh* pMesh, T& min, T& max, unsigned int dof=0);
        
    };


    template<typename T,typename I>
    DVector<T,I>::DVector()
    {
        m_uiData = NULL;
        m_uiSize =0;
        m_uiDof=0;
        m_uiIsElemental = false;
        m_uiIsGhosted =false;
        m_uiIsUnzip =false;
        m_uiComm = MPI_COMM_NULL;

    }

    template<typename T,typename I>
    DVector<T,I>::~DVector()
    {
        // delete [] m_uiData;
        // m_uiData = NULL;
        // m_uiSize =0;
        // m_uiDof=0;
        // m_uiIsElemental = false;
        // m_uiIsGhosted =false;
        // m_uiIsUnzip =false;
        // m_uiComm = MPI_COMM_NULL;
    }

    template<typename T,typename I>
    void DVector<T,I>::VecCreate(const ot::Mesh* pMesh, bool isGhosted, bool isUnzip, bool isElemental, unsigned int dof)
    {
        
        m_uiIsGhosted = isGhosted;
        m_uiIsUnzip = isUnzip;
        m_uiIsElemental = isElemental;
        m_uiDof = dof;
        m_uiComm = pMesh->getMPICommunicator();

        if(!(pMesh->isActive()))
        {
            m_uiData = NULL;
            return;
        }

        

        if(m_uiIsUnzip)
        {  
            //todo: make the unzip vector sizes correct for local and ghosted element vec. 
            if(m_uiIsGhosted)
            {
                if(m_uiIsElemental)
                    m_uiSize = pMesh->getLocalBlockList().size()*m_uiDof;
                else
                    m_uiSize = pMesh->getDegOfFreedomUnZip()*m_uiDof;
            }else
            {
                if(m_uiIsElemental)
                    m_uiSize = pMesh->getLocalBlockList().size()*m_uiDof;
                else
                    m_uiSize = pMesh->getDegOfFreedomUnZip()*m_uiDof;


            }
                

        }else
        {
            if(m_uiIsGhosted)
            {
                if(m_uiIsElemental)
                    m_uiSize = pMesh->getAllElements().size()*m_uiDof;
                else
                    m_uiSize = pMesh->getDegOfFreedom()*m_uiDof;
            }else
            {
                if(m_uiIsElemental)
                    m_uiSize = pMesh->getNumLocalMeshElements()*m_uiDof;
                else
                    m_uiSize = pMesh->getNumLocalMeshNodes()*m_uiDof;


            }

        }

        m_uiData = new T[m_uiSize];

        return ;

    }

    template<typename T, typename I>
    void DVector<T,I>::VecCreateDG(const ot::Mesh* pMesh, bool isGhosted, unsigned int dof)
    {
        m_uiIsGhosted = isGhosted;
        m_uiIsElemental = true; // dg vectors are inherently elemental vectors. 
        m_uiDof = dof;
        m_uiComm = pMesh->getMPICommunicator();

        if(!(pMesh->isActive()))
        {
            m_uiData = NULL;
            return;
        }

        if(m_uiIsGhosted)
            m_uiSize = pMesh->getAllElements().size()*pMesh->getNumNodesPerElement() * m_uiDof;
        else
            m_uiSize = pMesh->getNumLocalMeshElements()*pMesh->getNumNodesPerElement() * m_uiDof;

        m_uiData = new T[m_uiSize];
        return ;


    }


    template<typename T,typename I>
    void DVector<T,I>::VecDestroy()
    {
        delete [] m_uiData;
        m_uiData = NULL;
        m_uiSize =0;
        m_uiDof=0;
        m_uiIsElemental = false;
        m_uiIsGhosted =false;
        m_uiIsUnzip =false;
    }


    template<typename T, typename I>
    void DVector<T,I>::operator=(DVector<T,I> other)
    {
        m_uiData = other.GetVecArray();;
        m_uiDof = other.GetDof();
        m_uiSize = other.GetSize();
        m_uiIsElemental = other.IsElemental();
        m_uiIsGhosted = other.IsGhosted();
        m_uiIsUnzip = other.IsUnzip();
        m_uiComm = other.GetMPIComm();

        return;
    }


    template<typename T, typename I>
    void DVector<T,I>::Get2DArray(T**& v2d, bool isAlloc)
    {

        assert((m_uiSize%m_uiDof)==0);
        const I sz_per_dof = m_uiSize/m_uiDof;
        if(!isAlloc)
        {
            v2d = new T*[m_uiDof];
        }

        for(unsigned int i=0; i< m_uiDof; i++)
            v2d[i] = m_uiData + i*sz_per_dof;

        return;

    }


    template<typename T, typename I>
    void DVector<T,I>::Get2DVec(T** v2d)
    {

        assert((m_uiSize%m_uiDof)==0);
        const I sz_per_dof = m_uiSize/m_uiDof;
        
        for(unsigned int i=0; i< m_uiDof; i++)
            v2d[i] = m_uiData + i*sz_per_dof;

        return;

    }


    template<typename T, typename I >
    void DVector<T,I>::VecCopy(DVector<T,I> v, bool isAlloc)
    {
        if(!isAlloc)
        {
            m_uiData = new T [ v.GetSize() ];
        }else
        {
            assert(m_uiSize ==v.GetSize());
        }

        m_uiSize = v.GetSize();
        m_uiDof = v.GetDof();
        m_uiIsElemental = v.IsElemental();
        m_uiIsGhosted = v.IsGhosted();
        m_uiIsUnzip = v.IsUnzip();

        T* dptr = v.GetVecArray();
         
        std::memcpy(m_uiData, dptr, sizeof(T)*m_uiSize);
        return;
            
    }

    template<typename T, typename I>
    void DVector<T,I>::VecFMA(const ot::Mesh * pMesh , DVector<T,I> vec, T a, T b, bool localOnly)
    {   
        
        if(!(pMesh->isActive()))
            return;

        assert((this->IsElemental() == vec.IsElemental()) && (this->IsUnzip() == vec.IsUnzip()) && (this->IsGhosted() == vec.IsGhosted()) && (this->GetSize() == vec.GetSize()) && (this->GetDof() == vec.GetDof()) );

        const T* dptr = vec.GetVecArray();
        const I nPDOF = vec.GetSizePerDof();

        if(m_uiIsUnzip==true)
        {
            const ot::Block* const  blkList = pMesh->getLocalBlockList().data();
            const unsigned int numBlks = pMesh->getLocalBlockList().size();
            
            if(localOnly)
            {  
                // do the FMA only for block internal. 
                for (unsigned int v=0; v < m_uiDof; v++)
                {
                    for(unsigned int blk =0; blk< numBlks; blk ++)
                    {
                        const unsigned int offset = blkList[blk].getOffset();
                        const unsigned int nx = blkList[blk].getAllocationSzX();
                        const unsigned int ny = blkList[blk].getAllocationSzY();
                        const unsigned int nz = blkList[blk].getAllocationSzZ();
                        const unsigned int pw = blkList[blk].get1DPadWidth();

                        for(unsigned int k =pw; k < (nz-pw); k++)
                         for(unsigned int j=pw; j < (ny-pw); j++)
                          for(unsigned int i=pw; i< (nx-pw); i++)
                            m_uiData[v*nPDOF + offset + k*ny*nx + j* nx + i] = a*m_uiData[v*nPDOF + offset + k*ny*nx + j* nx + i] + b*dptr[v*nPDOF + offset + k*ny*nx + j* nx + i];
                        
                    }

                }


            }else
            {
                // do FMA operation to padding points as well. 
                for(I i=0; i< m_uiSize; i++)
                    m_uiData[i] = a*m_uiData[i] + b*dptr[i];

                
            }



        }else
        {
            if(localOnly)
            {
                for (unsigned int v=0; v < m_uiDof; v++)
                {
                    for(unsigned int n = pMesh->getNodeLocalBegin(); n < pMesh->getNodeLocalEnd(); n++ )
                        m_uiData[ v*nPDOF + n ] = a*m_uiData[v*nPDOF + n] + b* dptr[v*nPDOF + n];
                    
                }   
    
            }else
            {
                for (unsigned int v=0; v < m_uiDof; v++)
                {
                    for(unsigned int n = pMesh->getNodePreGhostBegin(); n < pMesh->getNodePostGhostEnd(); n++ )
                        m_uiData[ v*nPDOF +  n] = a*m_uiData[v*nPDOF +  n] + b* dptr[v*nPDOF +  n];

                }
                

            }

        }
        
        


    }

    template<typename T, typename I>
    void DVector<T,I>::VecFMA(const ot::Mesh * pMesh , T a, T b, bool localOnly)
    {
        if(!(pMesh->isActive()))
            return;


        const I nPDOF = this->GetSizePerDof();

        if(m_uiIsUnzip == true)
        {

            const ot::Block* const  blkList = pMesh->getLocalBlockList().data();
            const unsigned int numBlks = pMesh->getLocalBlockList().size();
            
            if(localOnly)
            {  
                // do the FMA only for block internal. 
                for (unsigned int v=0; v < m_uiDof; v++)
                {
                    for(unsigned int blk =0; blk< numBlks; blk ++)
                    {
                        const unsigned int offset = blkList[blk].getOffset();
                        const unsigned int nx = blkList[blk].getAllocationSzX();
                        const unsigned int ny = blkList[blk].getAllocationSzY();
                        const unsigned int nz = blkList[blk].getAllocationSzZ();
                        const unsigned int pw = blkList[blk].get1DPadWidth();

                        for(unsigned int k =pw; k < (nz-pw); k++)
                         for(unsigned int j=pw; j < (ny-pw); j++)
                          for(unsigned int i=pw; i< (nx-pw); i++)
                            m_uiData[v*nPDOF + offset + k*ny*nx + j* nx + i] = a * m_uiData[v*nPDOF + offset + k*ny*nx + j* nx + i] + b;
                        
                    }

                }


            }else
            {
                // do FMA operation to padding points as well. 
                for(I i=0; i< m_uiSize; i++)
                    m_uiData[i] = a * m_uiData[i] + b;

            }

        }else
        {
            if(localOnly)
            {   
                for(unsigned int v =0; v < m_uiDof; v++)
                {
                    for(unsigned int n = pMesh->getNodeLocalBegin(); n < pMesh->getNodeLocalEnd(); n++ )
                        m_uiData[ v*nPDOF +   n] = a*m_uiData[v*nPDOF +   n] + b;
                }
                
            }else
            {
                for(unsigned int v =0; v < m_uiDof; v++)
                {
                    for(unsigned int n = pMesh->getNodePreGhostBegin(); n < pMesh->getNodePostGhostEnd(); n++ )
                        m_uiData[v*nPDOF + n] = a*m_uiData[ v*nPDOF + n] + b;

                }
                
            }

        }
        
        
    }

    template<typename T, typename I>
    void DVector<T,I>::VecMinMax(const ot::Mesh* pMesh, T& min, T& max, unsigned int dof)
    {
        if(!(pMesh->isActive())) 
        {
            min=0; max=0;
            return;
        }

        assert( (m_uiIsUnzip==false) && (m_uiIsGhosted==true) && (m_uiIsElemental ==false) ); // current min max implementation. 
        const I sz = pMesh->getDegOfFreedom();
        T* ptr = m_uiData + dof*sz;

        const unsigned int num_local_nodes = pMesh->getNumLocalMeshNodes();
        MPI_Comm comm_active = pMesh->getMPICommunicator();
        min = vecMin(ptr + pMesh->getNodeLocalBegin(), num_local_nodes,comm_active);
        max = vecMax(ptr + pMesh->getNodeLocalBegin(), num_local_nodes,comm_active);

        return;

    }





}// end of namespace ot


typedef ot::DVector<DendroScalar, DendroIntL> DVec;