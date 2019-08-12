/**
 * @file nuts.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Spatially adaptive non uniform time stepping framework. 
 * @version 0.1
 * @date 2019-07-12
 * @copyright Copyright (c) 2019, School of Computing, University of Utah. 
 * 
 */
#pragma once
#include "mesh.h"
#include <functional>
#include "block.h"
#include "point.h"
#include <assert.h>


namespace ts
{   /**
    * @brief Basic: class for performing non-uniform time stepping. 
    * In order to perform non uniform time stepping the octree to block decomposition
    * needed to be completed. We assume that the blocks are setup in the mesh class. 
    */
    class NUTS
    {
        private:

            /**@brief: pointer to the octree mesh */
            ot::Mesh*  m_uiOctGrid;

            /**@brief: pointer to the block mesh */
            ot::Mesh* m_uiBlockGrid;

            /**@brief: number of time  integrator stages */
            unsigned int m_uiS;

            /**@brief: starting time value for simulation*/
            double m_uiTB;

            /**@brief: Ending time value for simulation*/
            double m_uiTE;

            /**@brief: dt minimum based on level */
            double m_uiDTMin;

            /**@brief: dt maximum based on level*/
            double m_uiDTMax;

            /**@brief: CFL factor  */
            double m_uiCFL;

            /**@brief: MPI active comm. */
            MPI_Comm m_uiCommActive;

            /**@brief: MPI global comm. */
            MPI_Comm m_uiCommGlobal;

            /**@brief: true if the comm. is active. */
            bool isActive;

            /**@brief: Domain point min*/
            Point m_uiPtMin;

            /**@brief: Domain point max*/
            Point m_uiPtMax;
            
            /**@brief: current time of the evolution */
            double m_uiCT;

            /**@brief: send scatter map level tag needed for sync operations */
            std::vector<unsigned int> m_uiSendSMLevelTag;

            /**@brief: recv scatter map level tag needed for sync operations */
            std::vector<unsigned int> m_uiRecvSMLevelTag;

            /**@brief: send counts for each level sync operation */
            unsigned int ** m_uiSendCountByLev;
            
            /**@brief: recv counts for each level sync operation */
            unsigned int ** m_uiRecvCountByLev;

            /**@brief: send offsets for each level sync operation */
            unsigned int ** m_uiSendOffsetByLev;

            /**@brief: recv offsets for each level sync operation */
            unsigned int ** m_uiRecvOffsetByLev;

            /**@brief: send reverse scatter map. */
            std::vector< std::pair<unsigned int,unsigned int> > m_uiSendNodeRSM;

            /**@brief: recv reverse scatter map. */
            std::vector< std::pair<unsigned int,unsigned int> > m_uiRecvNodeRSM;

            

            

        private:
            /** @brief perform required flag operations required for non-uniform timesteping */
            void compBlockFlags();

            /** @brief: Computes the block wise scatter map */
            void computeBlockSM();




        public:

            /**@brief: constructor
             * @param[in] octMesh: octree mesh
             * @param[in] s: number of rk stages. 
             * @param[in] tb: time begin
             * @param[in] te: time end
             * @param[in] cfl: cfl factor
             */
            NUTS(ot::Mesh* octMesh, unsigned int s, double tb, double te, Point ptMin, Point ptMax, double cfl=0.25);

            /**@brief: default destructor */
            ~NUTS();

            /**@brief: converts an zip vector to unzip format. 
             * @param[in] in: input zip vector
             * @param[out] out : output unzip vector
             * @param[in] isGSync: true if the ghost sync is performed
             * @param[in] dof: number of degrees of freedoms. 
             */
            template<typename T>
            void zip_to_unzip(T** in, T** out, bool isGSync=false, unsigned int dof=1);

            /** 
             * @brief: converts the unzip vector to zip vector
             * @param[in] in: input unzip vector 
             * @param[out] out: output zip vector
             * @param[in] isGsync : perform ghost sync after zip operation
             * @param[in] dof: number of degrees of freedoms
            */
            template<typename T>
            void unzip_to_zip(T** in, T** out, bool isGSync=false, unsigned int dof=1);


            /**
             * @brief: We need to perform block level sync operation
             * @param[in] in: input vector of multiple vars (unzip version). 
             * @param[in] blk: block ID, (should be sync on local blocks only. )
             * @param[in] dof: number of variables
             *  
            */
            template<typename T>
            void unzip_sync(T** in, unsigned int blk , unsigned int dof=1);
            
            /**
             * @brief: Sync operation accross all the blocks.  
             * @param [in] in: input vector
             * @param [in] dof: number of vars
            */
            template<typename T>
            void unzip_sync_all(T**in, unsigned int dof=1);

            /**
             * @brief: sync to the next coarsest time step  
            */
            void sync_ts();

            /**
             * @brief: perform single time step 
             * @param[in] in : input vector (previous time step)
             * @param[out] out: output vector (next timestep)
             * @param[in] f_rhs: functional to compute the right hand side of the PDE(s)
             * @param[in] a: time integrator coefficient matrix
             * @param[in] b: different stage combination coefficients
             * @param[in] c: different time stage coefficients
             * @param[in] dof: number of degrees of freedoms
             * 
            */
            template<typename T>
            void timestep(const T** in, T** out, std::function<void(T,T**,T**)>f_rhs, const double*a, const double *b, const double* c, unsigned int dof=1);            



        
    };

    template<typename T>
    void NUTS::zip_to_unzip(T** in, T** out, bool isGSync, unsigned int dof)
    {

        if(!isGSync)
        { 
          for(unsigned int v=0;v<dof;v++)
            m_uiOctGrid->performGhostExchange(in[v]);

        } 
        
        for(unsigned int v=0;v<dof;v++)
          m_uiOctGrid->unzip(in[v],out[v]);

    }


    template<typename T>
    void NUTS::unzip_to_zip(T** in, T** out, bool isGSync, unsigned int dof)
    {
        for(unsigned int v=0;v<dof;v++)
          m_uiOctGrid->zip(in[v],out[v]);

        if(isGSync)
        {
            for(unsigned int v=0;v<dof;v++)
                m_uiOctGrid->performGhostExchange(in[v]);
        }
    }

    template<typename T>
    void NUTS::unzip_sync(T** in, unsigned int blk , unsigned int dof)
    {
        
        const unsigned int blockLocalBegin = m_uiBlockGrid ->getElementLocalBegin();
        const unsigned int blockLocalEnd   = m_uiBlockGrid->getElementLocalEnd();

        if(blk< blockLocalBegin  || blk>= blockLocalEnd)
            return;
        
        //const unsigned int ;


    }

    template<typename T>
    void NUTS::unzip_sync_all(T**in, unsigned int dof)
    {
        const unsigned int blockLocalBegin = m_uiBlockGrid ->getElementLocalBegin();
        const unsigned int blockLocalEnd   = m_uiBlockGrid->getElementLocalEnd();

        for(unsigned int blk = blockLocalBegin; blk<blockLocalEnd ; blk++)
            unzip_sync(in,blk,dof);

    }



}


