/**
 * @file enuts.h
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
#include "ts.h"
#include "ets.h"
#include "meshUtils.h"
#include <bitset>
#include "subSM.h"
#include "oct2vtk.h"
#include <iostream>
#include "blkAsync.h"

namespace ts
{   


    /**
     * @brief : simple structure to support storing of a single time step (explicit methods)
     * note that the stages are numbered from 1 to m_uiNumStages. 
     * @tparam T : vector type
     */
    template <typename T>
    struct BlockTimeStep
    {

        public:
            /**@brief: ets stages
             * stage 0    : input vector
             * stage 1    : ets stage 1
             *            .
             *            .
             * stage p    : ets stage p
             * stage p+1  : output vector
            */
            std::vector< BlockAsyncVector<T> > _vec;

            /**@brief: rk stage*/
            unsigned int  _rks  = ETS_STAGE_DEFAULT;

            /**@brief: block time*/
            unsigned int  _time = LOOK_UP_TABLE_DEFAULT;

            /**
             * @brief allocate the Block time step vector. 
             * @param numVec : number of vectors per block. 
             * @param blkid  : block id
             * @param sz     : sizes of the each dimension
             * @param dof    : degrees of freedoms. 
             */
            void alloc_vec(unsigned int numVec, unsigned int blkid ,const unsigned int *sz , unsigned int dof=1)
            {
                _vec.resize(numVec);
                for(unsigned int k=0; k < _vec.size(); k++)
                    _vec[k].createVec(blkid,sz, false, BLK_ASYNC_VEC_MODE::BLK_UNZIP, dof);
            }

            /**
             * @brief deallocate the vectors. 
             */
            void dealloc_vec()
            {
                for(unsigned int k=0; k < _vec.size(); k++)
                    _vec[k].destroyVec();

                _vec.clear();
            }


    };


    
    /**
    * @brief Basic: class for performing non-uniform time stepping. (explicit time stepping)
    * In order to perform non uniform time stepping the octree to block decomposition
    * needed to be completed. We assume that the blocks are setup in the mesh class. 
    * 
    * Note: The stages are numbered from 1 to m_uiNumStages. 
    */
    #ifdef __PROFILE_ENUTS__
        enum ENUTSPROFILE {ENUTS_EVOLVE=0, BLK_SYNC, NUTS_CORRECTION, ENUTS_BLK_UNZIP, ENUTS_BLK_ZIP,  ENUTS_LAST};
    #endif

    template<typename T, typename Ctx>
    class ExplicitNUTS  : public ETS<T,Ctx>
    {

        
        using ETS<T,Ctx>::m_uiAppCtx;
        using ETS<T,Ctx>::m_uiNumStages;
        using ETS<T,Ctx>::m_uiEVar;
        using ETS<T,Ctx>::m_uiStVec;
        using ETS<T,Ctx>::m_uiEVecTmp;
        using ETS<T,Ctx>::m_uiTimeInfo;
        using ETS<T,Ctx>::m_uiAij;
        using ETS<T,Ctx>::m_uiBi;
        using ETS<T,Ctx>::m_uiCi;
        using ETS<T,Ctx>::m_uiType;
        using ETS<T,Ctx>::m_uiIsInternalAlloc;

        #ifdef __PROFILE_ENUTS__

            public:
                std::vector<profiler_t> m_uiPt = std::vector<profiler_t>(static_cast<int>(ENUTSPROFILE::ENUTS_LAST)); 
                const char *ENUTSPROFILE_NAMES[static_cast<int>(ENUTSPROFILE::ENUTS_LAST)] = {"evolve","blk_sync", "nuts_corr", "blk_unzip", "blk_zip" };

                void init_pt()
                {
                    for(unsigned int i=0; i < m_uiPt.size(); i++)
                        m_uiPt[i].start();

                    m_uiAppCtx->init_pt();
                }

                void reset_pt()
                {
                    for(unsigned int i=0; i < m_uiPt.size(); i++)
                        m_uiPt[i].snapreset();

                    m_uiAppCtx->reset_pt();
                }

                void dump_pt(std::ostream& outfile)
                {
                    const ot::Mesh* m_uiMesh = m_uiAppCtx->get_mesh();
                    
                    if(!(m_uiMesh->isActive()))
                        return;

                    int rank = m_uiMesh->getMPIRank();
                    int npes = m_uiMesh->getMPICommSize();

                    MPI_Comm comm = m_uiMesh->getMPICommunicator();
                    const unsigned int currentStep = m_uiAppCtx->get_ts_info()._m_uiStep;

                    double t_stat;
                    double t_stat_g[3];

                    if(!rank)
                    {
                        //writes the header
                        if(currentStep<=1)
                        outfile<<"step_nuts\t act_npes\t glb_npes\t maxdepth\t numOcts\t dof_cg\t dof_uz\t"<<\
                        "gele_min\t gele_mean\t gele_max\t"\
                        "lele_min\t lele_mean\t lele_max\t"\
                        "lnodes_min\t lnodes_mean\t lnodes_max\t"\
                        "remsh_igt_min\t remesh_igt_mean\t remesh_igt_max\t"\
                        "evolve_min\t evolve_mean\t evolve_max\t"\
                        "blk_sync_min\t blk_sync_mean\t blk_sync_max\t"\
                        "nuts_corr_min\t nuts_corr_mean\t nuts_corr_max\t"\
                        "blk_unzip_min\t blk_unzip_mean\t blk_unzip_max\t"\
                        "rhs_blk_min\t rhs_blk_mean\t rhs_blk_max\t"\
                        "blk_zip_min\t blk_zip_mean\t blk_zip_max\t"<<std::endl;
                    }

                        
                    if(!rank) outfile<<currentStep<<"\t ";
                    if(!rank) outfile<<m_uiMesh->getMPICommSize()<<"\t ";
                    if(!rank) outfile<<m_uiMesh->getMPICommSizeGlobal()<<"\t ";
                    if(!rank) outfile<<m_uiMaxDepth<<"\t ";


                    DendroIntL localSz=m_uiMesh->getNumLocalMeshElements();
                    DendroIntL globalSz;

                    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
                    if(!rank)outfile<<globalSz<<"\t ";

                    localSz=m_uiMesh->getNumLocalMeshNodes();
                    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
                    if(!rank)outfile<<globalSz<<"\t ";

                    localSz=m_uiMesh->getDegOfFreedomUnZip();
                    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
                    if(!rank)outfile<<globalSz<<"\t ";

                    DendroIntL ghostElements=m_uiMesh->getNumPreGhostElements()+m_uiMesh->getNumPostGhostElements();
                    DendroIntL localElements=m_uiMesh->getNumLocalMeshElements();

                    t_stat=ghostElements;
                    min_mean_max(&t_stat,t_stat_g,comm);
                    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

                    t_stat=localElements;
                    min_mean_max(&t_stat,t_stat_g,comm);
                    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

                    DendroIntL ghostNodes=m_uiMesh->getNumPreMeshNodes()+m_uiMesh->getNumPostMeshNodes();
                    DendroIntL localNodes=m_uiMesh->getNumLocalMeshNodes();

                    t_stat=localNodes;
                    min_mean_max(&t_stat,t_stat_g,comm);
                    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

                    t_stat= m_uiAppCtx->m_uiCtxpt[CTXPROFILE::REMESH].snap;
                    min_mean_max(&t_stat,t_stat_g,comm);
                    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

                    t_stat=m_uiPt[ENUTSPROFILE::ENUTS_EVOLVE].snap;
                    min_mean_max(&t_stat,t_stat_g,comm);
                    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

                    t_stat=m_uiPt[ENUTSPROFILE::BLK_SYNC].snap;
                    min_mean_max(&t_stat,t_stat_g,comm);
                    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

                    t_stat=m_uiPt[ENUTSPROFILE::NUTS_CORRECTION].snap;
                    min_mean_max(&t_stat,t_stat_g,comm);
                    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

                    t_stat=m_uiPt[ENUTSPROFILE::ENUTS_BLK_UNZIP].snap;
                    min_mean_max(&t_stat,t_stat_g,comm);
                    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

                    t_stat=m_uiAppCtx->m_uiCtxpt[CTXPROFILE::RHS_BLK].snap;
                    min_mean_max(&t_stat,t_stat_g,comm);
                    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

                    t_stat=m_uiPt[ENUTSPROFILE::ENUTS_BLK_ZIP].snap;
                    min_mean_max(&t_stat,t_stat_g,comm);
                    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

                    if(!rank) outfile<<std::endl;


                    





                }

        #endif
        
        protected:

            /**@brief: minimum level of the grid. */
            unsigned int m_uiLevMin;

            /**@brief: max level of the grid. */
            unsigned int m_uiLevMax;

            /**@brief: element to block map. */
            std::vector< unsigned int > m_uiE2B; 
            
            /**@brief: unzip vector for m_uiEVar*/
            DVec m_uiEvarUzip;

            /**@brief: DG vector for m_uiEVar*/
            DVec m_uiEvarDG;

            /**@brief: List of block async vector. */
            std::vector<ts::BlockTimeStep<T>> m_uiBVec;

            /**@brief : keep track of the element time. */
            std::vector<unsigned int> m_uiEleTime;
            
            /**@brief min and max time of blocks for a given level. */
            std::vector<unsigned int> m_uiBlkTimeLevMinMax;

            /**@brief: explicit timer interpolation operators for ENUTS */
            ENUTSOp* m_uiECOp = NULL;

            /**@brief: store the partial step id for the latest evolution. */
            unsigned int m_uiPt;

            

            


        private:

            /**@brief: allocates the data stuctures and initialize with the current mesh stores in the class. */
            void init_data_structures();

            /**@brief: freee the allocated data structures. */
            void free_data_structures();

            
        private:
            /**
             * @brief Allocates internal variables for the time stepper. 
             * @return int 
             */
            virtual int allocate_internal_vars();

            /**@brief: Deallocate internal variables*/
            virtual int deallocate_internal_vars(); 

            /**@bried: synchronize the block vectors*/
            virtual void sync_blk_timestep(unsigned int blk, unsigned int rk_s);
            
            /**@brief: update the element timestep*/
            virtual void update_ele_timestep();


            

        public:

            /**@brief: constructor
             * Assumptions: Note that blocks result in from octree to block decomposition can be not 2:1 balanced.
             * But to perform NUTS, blocks should be 2:1 balanced. 
            */
            ExplicitNUTS(Ctx* ctx);

            /**@brief: default destructor */
            ~ExplicitNUTS();

            /**@brief: build all the required data structures for the non uniform TS*/
            void init();

            /**@brief : evolve the variables to the next coarsest time step (i.e. loop over the block sequence. ) Evolution in the sense of the  NUTS*/
            virtual void evolve();

            /**
             * @brief : Perform parital step for the evolution. 
             * @param[in] pt: partial step id. 
             */
            void partial_evolve(unsigned int pt);

            /**@brief: compute the evolution with remesh triggered at partial steps. */
            void evolve_with_remesh(double sync_time, unsigned int remesh_freq=5);
            
            /**@brief: returns  a constant pointer to the sub scatter maps. */
            //inline ot::SubScatterMap* const get_sub_scatter_maps() const {return m_uiSubSM;}

            /**@brief: perform remesh for the nuts class. */
            virtual void remesh(unsigned int grain_sz = DENDRO_DEFAULT_GRAIN_SZ, double ld_tol = DENDRO_DEFAULT_LB_TOL, unsigned int sf_k = DENDRO_DEFAULT_SF_K);

            /**@brief: update all the internal data strucutures, with a new mesh pointer. */
            virtual int sync_with_mesh();

            /**@brief returns the dt min value */
            T get_dt_min() const { return m_uiTimeInfo._m_uiTh; } 
            
            /**@brief returns the dt max value */
            T get_dt_max() const { return m_uiTimeInfo._m_uiTh*(1u<<(m_uiLevMax-m_uiLevMin)); }

            /**@brief: prints the load balance statistis*/
            void  dump_load_statistics(std::ostream & sout) const ;

            /**@brief: computes the estimated speedup for a given mesh. */
            void  dump_est_speedup(std::ostream & sout);

            /**@brief : perform blk time step vector to vec CG zip operation. 
             * @param blkVec : pointer to the block vector. 
             * @param vecCG  : allocated CG vector. 
             * @param s: blkVec.__vec[s] access parameter. 
            */
            void blk_vec_to_zipCG(const ts::BlockTimeStep<T>* blkVec, DVec& vecCG,unsigned int s=0);

            /**
             * @brief perform zip operation (DG) to get the vecDG from block vector.
             * @param blkVec : blk vector in
             * @param vecDG  : allocated DG vector. 
             * @param s: blkVec.__vec[s] access parameter. 
             */
            void blk_vec_to_zipDG(const ts::BlockTimeStep<T>* blkVec, DVec& vecDG,unsigned int s=0);

            /**
             * @brief copy DG vector to the block vector. 
             * 
             * @param vecDG 
             * @param blkVec 
             * @param s 
             */
            void vecDG_to_blk_vec(DVec& vecDG, ts::BlockTimeStep<T>* blkVec, unsigned int s=0);
        
    };

    template<typename T, typename Ctx>
    ExplicitNUTS<T,Ctx>::ExplicitNUTS(Ctx* ctx) : ETS<T,Ctx>(ctx)
    {
        this->init_data_structures();
        
    }

    template<typename T, typename Ctx>
    ExplicitNUTS<T,Ctx>::~ExplicitNUTS()
    {
        this->deallocate_internal_vars();
        this->free_data_structures();

    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::dump_est_speedup(std::ostream & sout)  
    {
        const ot::Mesh* pMesh =  m_uiAppCtx->get_mesh();

        if(pMesh->isActive())
        {
            const unsigned int rank = pMesh->getMPIRank();
            const unsigned int npes = pMesh->getMPICommSize();
            MPI_Comm comm = pMesh->getMPICommunicator();

            const unsigned int  finest_t = 1u<<(m_uiLevMax - m_uiLevMax); // finest  time level.
            const unsigned int coarset_t = 1u<<(m_uiLevMax - m_uiLevMin); // coarset time level. 

            double red_stat[3]; // for mpi reduction stats. 
            std::vector<ot::Block> blkList = pMesh->getLocalBlockList();
            const ot::TreeNode* pNodes = pMesh->getAllElements().data();

            // initialize the block async vectors
            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                m_uiBVec[blk]._time = 0;
                m_uiBVec[blk]._rks  = 0;
            }


            double localSz = 0;
            
            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                double nnx = (blkList[blk].getAllocationSzX() -2*blkList[blk].get1DPadWidth() );
                localSz += nnx*nnx*nnx;  //(blkList[blk].getLocalElementEnd() - blkList[blk].getLocalElementBegin());
            }

            double globalESz;
            par::Mpi_Reduce(&localSz,&globalESz,1,MPI_SUM,0, comm);

            double globalBlkSz;
            localSz = pMesh->getLocalBlockList().size();
            par::Mpi_Reduce(&localSz,&globalBlkSz,1,MPI_SUM,0, comm);

            double totalW=0;
            for(unsigned int pt=0; pt< coarset_t; pt ++)
            {
                double numBlkEvolved = 0; 

                for(unsigned int blk =0; blk < blkList.size(); blk++)
                {
                    const unsigned int BLK_T = m_uiBVec[blk]._time;
                    const unsigned int NN   =  m_uiBVec[blk]._vec[0].getSz();
                    const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
                    const unsigned int BLK_DT = 1u<<(m_uiLevMax - bLev);
                    if( pt% BLK_DT !=0 )
                        continue;

                    double nnx = (blkList[blk].getAllocationSzX() -2*blkList[blk].get1DPadWidth() );
                    //std::cout<<"pt: "<<pt<<" blk id : "<<blk<<" blklevel : "<<bLev<<" is now at btime : "<<BLK_T<<" dt: "<<BLK_DT<<std::endl;

                    numBlkEvolved += nnx*nnx*nnx;  //(blkList[blk].getLocalElementEnd() - blkList[blk].getLocalElementBegin());

                }

                par::Mpi_Reduce(&numBlkEvolved,red_stat,1,MPI_MIN,0,comm);
                par::Mpi_Reduce(&numBlkEvolved,red_stat+1,1,MPI_SUM,0,comm);
                par::Mpi_Reduce(&numBlkEvolved,red_stat+2,1,MPI_MAX,0,comm);

                totalW += red_stat[1];

                if(!rank)
                    std::cout<<" partial step : "<<pt<< " number of dof evolved (min, sum, max): ("<<red_stat[0]<<",\t "<<red_stat[1]<<",\t "<<red_stat[2]<<") out of  "<<globalESz<<std::endl;



                // compute the time step vector and increment time. 
                for(unsigned int blk =0; blk < blkList.size(); blk++)
                {

                    const unsigned int BLK_T = m_uiBVec[blk]._time;
                    const unsigned int NN   =  m_uiBVec[blk]._vec[0].getSz();
                    const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
                    const unsigned int BLK_DT = 1u<<(m_uiLevMax - bLev);
                    if( pt% BLK_DT !=0 )
                        continue;

                    m_uiBVec[blk]._time += 1u<<(m_uiLevMax -bLev); 
                    m_uiBVec[blk]._rks=0;
                    
                
                }

            }

            if(!rank)
            {   
                std::cout<<"nuts W : "<<totalW<<" uts: "<<globalESz*coarset_t<<std::endl;
                std::cout<<" ets. speedup : "<<(globalESz*coarset_t/(double)(totalW))<<std::endl;
            }
                


            
        }




    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::init_data_structures()
    {

         // identify dependent and independent blocks. 
        const ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        const bool isActive = pMesh->isActive();
        
        pMesh->computeMinMaxLevel(m_uiLevMin,m_uiLevMax);
        par::Mpi_Bcast(&m_uiLevMin,1,0,pMesh->getMPIGlobalCommunicator());
        par::Mpi_Bcast(&m_uiLevMax,1,0,pMesh->getMPIGlobalCommunicator());
        assert( (m_uiLevMin > 0 ) && (m_uiLevMax <= m_uiMaxDepth) );

    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::free_data_structures()
    {
        
        return ;

    }
    
    template<typename T, typename Ctx>
    int ExplicitNUTS<T,Ctx>::allocate_internal_vars()
    {

        assert(m_uiNumStages>0);
        const unsigned int DOF = m_uiEVar.GetDof();

        if(m_uiIsInternalAlloc)
            return 0; // no need to allocated again if the internal vars are allocated. 

        m_uiStVec.resize(m_uiNumStages);
        
        // for(unsigned int i=0; i < m_uiNumStages; i++)
        //     m_uiStVec[i].VecCreate(m_uiAppCtx->get_mesh(), false , true, false , m_uiEVar.GetDof());

        for(unsigned int i=0; i < m_uiNumStages; i++)
            m_uiStVec[i].VecCreateDG(m_uiAppCtx->get_mesh(), true, m_uiEVar.GetDof());

        m_uiEVecTmp[0].VecCreate(m_uiAppCtx->get_mesh(), m_uiEVar.IsGhosted() , m_uiEVar.IsUnzip(), m_uiEVar.IsElemental() , m_uiEVar.GetDof());
        m_uiEVecTmp[1].VecCreate(m_uiAppCtx->get_mesh(), m_uiEVar.IsGhosted() , m_uiEVar.IsUnzip(), m_uiEVar.IsElemental() , m_uiEVar.GetDof());

        m_uiEvarUzip.VecCreate(m_uiAppCtx->get_mesh(), false , true, false ,m_uiEVar.GetDof());
        m_uiEvarDG.VecCreateDG(m_uiAppCtx->get_mesh(), true, m_uiEVar.GetDof());

        const ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        if(pMesh->isActive())
        {
            m_uiEleTime.resize(pMesh->getAllElements().size(),0);
            std::vector<ot::Block> blkList = pMesh->getLocalBlockList();
            m_uiBVec.resize(blkList.size());

            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                const unsigned int sz[3] = { blkList[blk].getAllocationSzX(), blkList[blk].getAllocationSzY(), blkList[blk].getAllocationSzZ()};
                m_uiBVec[blk].alloc_vec(m_uiNumStages+2, blk, sz, DOF);
            }

        }

        m_uiECOp = new ENUTSOp(m_uiType);
        m_uiIsInternalAlloc=true;
        m_uiBlkTimeLevMinMax.resize(2*m_uiMaxDepth,0);
        return 0;

        
    }

    template<typename T, typename Ctx>
    int ExplicitNUTS<T,Ctx>::deallocate_internal_vars()
    {

        if(!m_uiIsInternalAlloc)
            return 0;

        for(unsigned int i=0; i < m_uiNumStages; i++)
        {
            this->m_uiStVec[i].VecDestroy();
        }

        this->m_uiStVec.clear();
        this->m_uiEVecTmp[0].VecDestroy();
        this->m_uiEVecTmp[1].VecDestroy();

        m_uiEvarUzip.VecDestroy();
        m_uiEvarDG.VecDestroy();


        for(unsigned int k=0; k < m_uiBVec.size(); k++)
        {
            for(unsigned int j=0; j < m_uiBVec[k]._vec.size(); j++)
                m_uiBVec[k]._vec[j].destroyVec();
            
            m_uiBVec[k]._vec.clear();
        }
        
        m_uiBVec.clear();

        delete m_uiECOp;
        m_uiIsInternalAlloc = false;
        return 0;

    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::blk_vec_to_zipCG(const ts::BlockTimeStep<T>* blkVec, DVec& vecCG, unsigned int s)
    {
        ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        const std::vector<ot::Block> blkList = pMesh->getLocalBlockList();
        const unsigned int DOF = vecCG.GetDof();
        
        for(unsigned int blk =0; blk < blkList.size(); blk++)
            blkVec[blk]._vec[s].zip(pMesh, vecCG.GetVecArray(), DOF);
        
        return;

    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::blk_vec_to_zipDG(const ts::BlockTimeStep<T>* blkVec, DVec& vecDG, unsigned int s)
    {
        ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        const std::vector<ot::Block> blkList = pMesh->getLocalBlockList();
        const unsigned int DOF = vecDG.GetDof();
        
        for(unsigned int blk =0; blk < blkList.size(); blk++)
            blkVec[blk]._vec[s].zipDG(pMesh, vecDG.GetVecArray(), DOF);
        
        return;
    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::vecDG_to_blk_vec(DVec& vecDG, ts::BlockTimeStep<T>* blkVec, unsigned int s)
    {
        ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        const std::vector<ot::Block> blkList = pMesh->getLocalBlockList();
        const unsigned int DOF = vecDG.GetDof();
        

        for(unsigned int blk =0; blk < blkList.size(); blk++)
        {
            blkVec[blk]._vec[s].copyFromVecDG(pMesh,vecDG.GetVecArray(),DOF);
            blkVec[blk]._vec[s].mark_unsynced();
        }
            



    }


    template<typename T, typename Ctx>
    int ExplicitNUTS<T,Ctx>::sync_with_mesh()
    {
        if(m_uiAppCtx -> is_ets_synced())
            return 0;

        this->deallocate_internal_vars();
        this->free_data_structures();

        this->init_data_structures();
        this->allocate_internal_vars();

        m_uiEVar = m_uiAppCtx->get_evolution_vars();
        m_uiAppCtx->set_ets_synced(true);

        return 0;


    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::init()
    {

        m_uiAppCtx->initialize();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        allocate_internal_vars();

    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::sync_blk_timestep(unsigned int blk, unsigned int rk_s)
    {

        #ifdef __PROFILE_ENUTS__
            m_uiPt[ENUTSPROFILE::BLK_SYNC].start();
        #endif

        ot::Mesh* pMesh = (ot::Mesh*)m_uiAppCtx->get_mesh();
        if((!(pMesh->isActive())) || m_uiBVec[blk]._vec[rk_s].isSynced() )
            return;

        

        const ot::TreeNode* pNodes   =   pMesh->getAllElements().data();
        const ot::Block* blkList     =   pMesh->getLocalBlockList().data();

        const unsigned int regLevel  =   blkList[blk].getRegularGridLev();
        const ot::TreeNode blkNode   =   blkList[blk].getBlockNode();
        const unsigned int PW        =   blkList[blk].get1DPadWidth();


        const unsigned int eOrder    =   pMesh->getElementOrder();
        const unsigned int nPe       =   pMesh->getNumNodesPerElement();
        
        const unsigned int * etVec   =   m_uiEleTime.data();
        MPI_Comm comm = pMesh->getMPICommunicator();
        
        assert(rk_s>=1 && rk_s <= m_uiNumStages);
        assert(PW>0);
        
        
        const unsigned int lx     =  blkList[blk].getAllocationSzX();
        const unsigned int ly     =  blkList[blk].getAllocationSzY();
        const unsigned int lz     =  blkList[blk].getAllocationSzZ();
        const unsigned int offset =  blkList[blk].getOffset(); 

        const unsigned int dgSz  =  pMesh->getAllElements().size() * pMesh->getNumNodesPerElement();
        const unsigned int cgSz  =  pMesh->getDegOfFreedom();
        const unsigned int unSz  =  pMesh->getDegOfFreedomUnZip();   

        const unsigned int* e2n  =  pMesh->getE2NMapping().data();
        const unsigned int* e2e  =  pMesh->getE2EMapping().data();


        const unsigned int bLev  =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
        const unsigned int bTime =  etVec[blkList[blk].getLocalElementBegin()];
        
        
        T dt;
        unsigned int tl=0;
        unsigned int lookUp;
        const unsigned int cSz[3] = { eOrder+1, eOrder+1, eOrder+1 };
        
        unsigned int fchild[4];
        unsigned int echild[2];

        const unsigned int dof = m_uiEVar.GetDof();
        std::vector<T> cVec;
        cVec.resize(dof*nPe*rk_s);

        T* cVin[rk_s*dof];  // correction vec. pointer in,
        T* cVout[rk_s*dof]; // correction vec. pointer out,

        T * dgWVec = m_uiEvarDG.GetVecArray();
        T * cgWVec = m_uiEVecTmp[0].GetVecArray();
        T * uzWVec = m_uiEvarUzip.GetVecArray();

        T* m_uiVec = (T*) m_uiBVec[blk]._vec[rk_s].data();

        T * dgStages[m_uiNumStages];
        for(unsigned int i=0; i < m_uiNumStages; i++)
            dgStages[i] = m_uiStVec[i].GetVecArray();


        std::vector<unsigned int> eid;
        computeBlockUnzipDepElements(pMesh, blk, eid);
        

        // apply time interpolation correction operators. 
        for(unsigned int m=0; m < eid.size(); m++)
        {
            const unsigned int elem = eid[m];
            tl = etVec[elem];
            
            if(pNodes[elem].getLevel() == bLev)
            {
                
                if( tl != bTime)
                {
                    std::cout<<" blk : "<<blk <<" bTime : "<<bTime <<" neigh. element (same lev) : "<<pNodes[elem]<<"  is at time : "<<tl<<"   time gap should be zero  "<<std::endl;
                    MPI_Abort(comm, 0);
                }

                #pragma unroll
                for(unsigned int v=0; v < dof; v++)
                    std::memcpy(dgWVec + v*dgSz + elem * nPe , dgStages[rk_s-1] + v*dgSz + elem * nPe, sizeof(T)*nPe );

            }else if(pNodes[elem].getLevel() < bLev)
            {
                const T dt_c             =  (1u<<(m_uiLevMax-pNodes[elem].getLevel()))*(m_uiTimeInfo._m_uiTh);
                const T dt_f             =  0.5*dt_c;

                if( tl == bTime)
                    dt=0;
                else
                {
                    
                    if( tl != (bTime  + (1u<<(m_uiLevMax - bLev)) )  )
                    {
                        std::cout<<" blk : "<<blk <<" bTime : "<<bTime <<" neigh. element (coaser) : "<<pNodes[elem]<<" is at time : "<<tl<<" time gap is "<<(1u<<(m_uiLevMax - bLev))<<" invalid" <<std::endl;
                        MPI_Abort(comm, 0);
                    }
                    dt=0.5*dt_c;
                    

                }

                for(unsigned int v=0; v < dof; v++)
                for(unsigned int s=1; s <= rk_s; s++ )
                {
                    cVin [v*rk_s  +  (s-1) ]  = dgStages[s-1] + v *dgSz  + elem*nPe;
                    cVout[v*rk_s  +  (s-1) ]  = cVec.data() + v*rk_s*nPe + (s-1)*nPe; 
                }

                #ifdef __PROFILE_ENUTS__
                    m_uiPt[ENUTSPROFILE::NUTS_CORRECTION].start();
                #endif

                m_uiECOp->Ccf(cVout,(const T**) cVin,cSz,rk_s,dt_c,dt_f,dt,dof);

                #ifdef __PROFILE_ENUTS__
                    m_uiPt[ENUTSPROFILE::NUTS_CORRECTION].stop();
                #endif
                
                #pragma unroll
                for(unsigned int v=0; v < dof; v++)
                    std::memcpy(dgWVec + v*dgSz + elem * nPe , cVout[v*rk_s  +  (rk_s-1)], sizeof(T)*nPe );

            }
            else
            {

                const T dt_c             =  (1u<<(m_uiLevMax-bLev))*(m_uiTimeInfo._m_uiTh);
                const T dt_f             =  0.5*dt_c;

                assert(pNodes[elem].getLevel() == bLev + 1);
                if( tl == bTime)
                    dt=0;
                else
                {
                    
                    if( bTime != (tl  + (1u<<(m_uiLevMax - bLev-1))) )
                    {
                        std::cout<<" blk : "<<blk <<" bTime : "<<bTime <<" neigh. element (finer) : "<<pNodes[elem]<<" is at time : "<<tl<<" time gap is "<<(1u<<(m_uiLevMax - bLev))<<" invalid" <<std::endl;
                        MPI_Abort(comm, 0);
                    }

                    dt=0.5*dt_c;
                    

                }

                for(unsigned int v=0; v < dof; v++)
                for(unsigned int s=1; s <= rk_s; s++ )
                {
                    cVin [v*rk_s  +  (s-1) ]  = dgStages[s-1] + v *dgSz  + elem*nPe;
                    cVout[v*rk_s  +  (s-1) ]  = cVec.data() + v*rk_s*nPe + (s-1)*nPe; 
                }

                #ifdef __PROFILE_ENUTS__
                    m_uiPt[ENUTSPROFILE::NUTS_CORRECTION].start();
                #endif
                
                m_uiECOp->Cfc(cVout,(const T**) cVin,cSz,rk_s,dt_c,dt_f,dt,dof);
                
                #ifdef __PROFILE_ENUTS__
                    m_uiPt[ENUTSPROFILE::NUTS_CORRECTION].stop();
                #endif
                
                #pragma unroll
                for(unsigned int v=0; v < dof; v++)
                    std::memcpy(dgWVec + v*dgSz + elem * nPe , cVout[v*rk_s  +  (rk_s-1)], sizeof(T)*nPe );

            }

        }

        // for(unsigned int elem = blkList[blk].getLocalElementBegin(); elem < blkList[blk].getLocalElementEnd(); elem++)
        // {
        //     const unsigned int ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLevel);
        //     const unsigned int ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLevel);
        //     const unsigned int ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLevel);

        //     const unsigned int emin = 0;
        //     const unsigned int emax = (1u<<(regLevel-blkNode.getLevel()))-1;

        //     assert(etVec[elem] == bTime );

        //     #pragma unroll
        //     for(unsigned int v=0; v < dof; v++)
        //         std::memcpy(dgWVec + v*dgSz + elem * nPe , dgStages[rk_s-1] + v*dgSz + elem * nPe, sizeof(T)*nPe );

        // }

        // for(unsigned int elem = blkList[blk].getLocalElementBegin(); elem < blkList[blk].getLocalElementEnd(); elem++)
        //     eid.push_back(elem);
        // pMesh->DG2CGVec(dgWVec,cgWVec,true,eid.data(),eid.size(),dof);
        // pMesh->unzip(cgWVec,uzWVec, &blk, 1, dof);
        
        // #pragma unroll
        // for(unsigned int v=0; v < dof; v++)
        //     std::memcpy(m_uiVec + v*lx*ly*lz , uzWVec + v* unSz + offset, sizeof(T)*(lx*ly*lz));

        // return;

        
        // now need to copy to the block unzip/ block asyncVector 
        const double  hx = (1u<<(m_uiMaxDepth-bLev))/(double)eOrder;
        
        const double xmin = blkNode.minX() - PW*hx; const double xmax = blkNode.maxX() + PW*hx;
        const double ymin = blkNode.minY() - PW*hx; const double ymax = blkNode.maxY() + PW*hx;
        const double zmin = blkNode.minZ() - PW*hx; const double zmax = blkNode.maxZ() + PW*hx;

        std::vector<ot::TreeNode> childOct;
        childOct.reserve(NUM_CHILDREN);

        std::vector<T> p2cI;
        p2cI.resize(nPe);

        for(unsigned int m=0; m < eid.size(); m++)
        {
            const unsigned int ele = eid[m];
            
            // no interpolation needed just copy. 
            if(pNodes[ele].getLevel()==bLev)
            {
                const double hh = (1u<<(m_uiMaxDepth - pNodes[ele].getLevel()))/(double) eOrder;
                const double invhh = 1.0/hh;
                
                for(unsigned int k=0; k < eOrder+1; k++)
                {
                    const double zz  = pNodes[ele].minZ() + k*hh;
                    if(zz < zmin || zz > zmax) 
                        continue;
                    const unsigned int kkz = std::round((zz-zmin)*invhh);
                    assert( std::fabs(zz-zmin-kkz*hh) < 1e-6);
                    assert(kkz >= 0 && kkz < lz);

                    for(unsigned int j=0; j < eOrder+1; j++)
                    {   
                        const double yy  = pNodes[ele].minY() + j*hh;
                        if(yy < ymin || yy > ymax) 
                            continue;
                        const unsigned int jjy = std::round((yy-ymin)*invhh);
                        //std::cout<<"yy: "<<yy<<" (ymin + hh*jjy): "<<(ymin + hh*jjy)<<std::endl;
                        assert( std::fabs(yy-ymin-jjy*hh) < 1e-6);
                        assert(jjy>=0 && jjy<ly);

                        for(unsigned int i=0; i < eOrder+1; i++)
                        {
                            const double xx = pNodes[ele].minX() + i*hh;
                            
                            if(xx < xmin || xx > xmax) 
                                continue;
                            const unsigned int iix = std::round((xx-xmin)*invhh);
                            assert( std::fabs(xx-xmin-iix*hh) < 1e-6);
                            assert(iix>=0 && iix<lx);

                            //std::cout<<"blk: "<<blk<<" copy : (i,j,k): ("<<kkz<<" , "<<jjy<<", "<<iix<<")"<<" of : "<<lx<<std::endl;

                            // if(blkNode.isAncestor(pNodes[ele]))
                            // {
                            //     const unsigned int ei=(pNodes[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLevel);
                            //     const unsigned int ej=(pNodes[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLevel);
                            //     const unsigned int ek=(pNodes[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLevel);

                            //     std::cout<<"blk: "<<blk<<" copy : (i,j,k): ("<<iix<<" , "<<jjy<<", "<<kkz<<")"<<" of : "<<lx<<" ei: "<<ei<<" ej: "<<ej<<" ek: "<<ek << " xx: "<<xx<<" yy: "<<yy<<" zz:"<<zz<<" xmin: "<<xmin<<" ymin: "<<ymin<<" zmin: "<<zmin<<" hh : "<<hh<<" hhx : "<<hx<<std::endl;

                            // }

                            for(unsigned int v=0; v < dof; v++)
                                m_uiVec[(v*lx*ly*lz) + kkz*lx*ly + jjy*lx + iix] =  dgWVec[v*dgSz + ele*nPe + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];

                        }
                    
                    }
                
                }


            }else if(pNodes[ele].getLevel() > bLev)
            {
                assert((bLev+1) == pNodes[ele].getLevel());
                const unsigned int cnum = pNodes[ele].getMortonIndex();
                ot::TreeNode tmpParent = pNodes[ele].getParent();
                
                const double hh = (1u<<(m_uiMaxDepth - pNodes[ele].getLevel()))/(double) eOrder;
                const double invhh = 1.0/(2*hh);

                assert(eOrder>1);
                const unsigned int cb =(eOrder%2==0) ? 0 : 1;

                for(unsigned int k=cb; k < eOrder+1; k+=2)
                {
                    const double zz  = (pNodes[ele].minZ() + k*hh);
                    if(zz < zmin || zz > zmax) 
                        continue;
                    const unsigned int kkz = std::round((zz-zmin)*invhh);
                    assert(kkz >= 0 && kkz < lz);

                    for(unsigned int j=cb; j < eOrder+1; j+=2)
                    {   
                        const double yy  = pNodes[ele].minY() + j*hh;
                        if(yy < ymin || yy > ymax) 
                            continue;

                        const unsigned int jjy = std::round((yy-ymin)*invhh);
                        assert(jjy>=0 && jjy<ly);

                        for(unsigned int i=cb; i < eOrder+1; i+=2)
                        {
                            const double xx = pNodes[ele].minX() + i*hh;
                            
                            if(xx < xmin || xx > xmax) 
                                continue;
                            const unsigned int iix = std::round((xx-xmin)*invhh);
                            assert(iix>=0 && iix<lx);


                            // if(blkNode.isAncestor(pNodes[ele]))
                            // {
                            //     const unsigned int ei=(pNodes[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLevel);
                            //     const unsigned int ej=(pNodes[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLevel);
                            //     const unsigned int ek=(pNodes[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLevel);

                            //     std::cout<<"blk: "<<blk<<" copy : (i,j,k): ("<<iix<<" , "<<jjy<<", "<<kkz<<")"<<" of : "<<lx<<" ei: "<<ei<<" ej: "<<ej<<" ek: "<<ek << " xx: "<<xx<<" yy: "<<yy<<" zz:"<<zz<<" xmin: "<<xmin<<" ymin: "<<ymin<<" zmin: "<<zmin<<" hh : "<<hh<<" hhx : "<<hx<<std::endl;

                            // }

                            //std::cout<<"blk: "<<blk<<" blk copy : (i,j,k): ("<<iix<<" , "<<jjy<<", "<<kkz<<")"<<" of : "<<lx<<" xx: "<<xx<<" yy: "<<yy<<" zz:"<<zz<<" xmin: "<<xmin<<" ymin: "<<ymin<<" zmin: "<<zmin<<" hh : "<<hh<<" hhx : "<<hx<<std::endl;
                            for(unsigned int v=0; v < dof; v++)
                                m_uiVec[(v*lx*ly*lz) + kkz*lx*ly + jjy*lx + iix] = dgWVec[v*dgSz + ele*nPe + k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];

                        }
                    
                    }
                
                }


              
            }else
            {
                assert((bLev) == (pNodes[ele].getLevel()+1));
                childOct.clear();
                pNodes[ele].addChildren(childOct);

                for(unsigned int child = 0; child < NUM_CHILDREN; child++)
                {

                    if( (childOct[child].maxX() < xmin || childOct[child].minX() >=xmax)  || (childOct[child].maxY() < ymin || childOct[child].minY() >=ymax) || (childOct[child].maxZ() < zmin || childOct[child].minZ() >=zmax) )
                        continue;


                    //std::cout<<"blk: "<<blk<<" blkNode: "<<blkNode<<" child: "<<child<<" child node "<<childOct[child]<<" parent : "<<pNodes[ele]<<std::endl;
                    const double hh = (1u<<(m_uiMaxDepth - childOct[child].getLevel()))/(double) eOrder;
                    const double invhh = 1.0/hh;

                    for(unsigned int v=0; v < dof; v++)
                    {
                        pMesh->parent2ChildInterpolation(&dgWVec[v*dgSz + ele*nPe],p2cI.data(),child,m_uiDim);

                        for(unsigned int k=0; k < eOrder+1; k++)
                        {
                            const double zz  = childOct[child].minZ() + k*hh;
                            if(zz < zmin || zz > zmax) 
                                continue;
                            const unsigned int kkz = std::round((zz-zmin)*invhh);
                            assert(kkz >= 0 && kkz < lz);

                            for(unsigned int j=0; j < eOrder+1; j++)
                            {   
                                const double yy  = childOct[child].minY() + j*hh;
                                if(yy < ymin || yy > ymax) 
                                    continue;

                                const unsigned int jjy = std::round((yy-ymin)*invhh);
                                assert(jjy>=0 && jjy<ly);

                                for(unsigned int i=0; i < eOrder+1; i++)
                                {
                                    const double xx = childOct[child].minX() + i*hh;
                                    
                                    if(xx < xmin || xx > xmax) 
                                        continue;
                                    const unsigned int iix = std::round((xx-xmin)*invhh);
                                    assert(iix>=0 && iix<lx);
                                    //std::cout<<"blk: "<<blk<<" blk copy : (i,j,k): ("<<iix<<" , "<<jjy<<", "<<kkz<<")"<<" of : "<<lx<<" xx: "<<xx<<" yy: "<<yy<<" zz:"<<zz<<" xmin: "<<xmin<<" ymin: "<<ymin<<" zmin: "<<zmin<<" hh : "<<hh<<" hhx : "<<hx<<" child: "<<child<<std::endl;
                                    m_uiVec[(v*lx*ly*lz) + kkz*lx*ly + jjy*lx + iix] =  p2cI[k*(eOrder+1)*(eOrder+1)+ j*(eOrder+1) + i];

                                }
                            
                            }
                        
                        }
                        
                    }

                }

            }



        }
        

        for(unsigned int elem = blkList[blk].getLocalElementBegin(); elem < blkList[blk].getLocalElementEnd(); elem++)
        {
            const unsigned int ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLevel);
            const unsigned int ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLevel);
            const unsigned int ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLevel);

            const unsigned int emin = 0;
            const unsigned int emax = (1u<<(regLevel-blkNode.getLevel()))-1;

            assert(etVec[elem] == bTime );

            // #pragma unroll
            // for(unsigned int v=0; v < dof; v++)
            //     std::memcpy(dgWVec + v*dgSz + elem * nPe , dgStages[rk_s-1] + v*dgSz + elem * nPe, sizeof(T)*nPe );
            
            for(unsigned int v=0; v < dof; v++)
             for(unsigned int k=0;k<(eOrder+1);k++)
              for(unsigned int j=0;j<(eOrder+1);j++)
               for(unsigned int i=0;i<(eOrder+1);i++)
                m_uiVec[v*lx*ly*lz + (ek*eOrder+k+PW)*(ly*lx)+(ej*eOrder+j+PW)*(lx)+(ei*eOrder+i+PW)]=dgStages[rk_s-1][ v*dgSz + elem * nPe + k*(eOrder+1)*(eOrder+1) + j*(eOrder+1) + i];

        }

        #ifdef __PROFILE_ENUTS__
            m_uiPt[ENUTSPROFILE::BLK_SYNC].stop();
        #endif


    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::update_ele_timestep()
    {
        ot::Mesh* pMesh = m_uiAppCtx->get_mesh();

        if(!(pMesh->isActive()))
            return;

        std::vector<ot::Block> blkList = pMesh->getLocalBlockList();
        for(unsigned int blk =0; blk < blkList.size(); blk++)
        {
            for(unsigned int ele = blkList[blk].getLocalElementBegin(); ele < blkList[blk].getLocalElementEnd(); ele ++)
                m_uiEleTime[ele] = m_uiBVec[blk]._time;
            
        }

    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::partial_evolve(unsigned int pt)
    {
        
        ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        
        if(!(pMesh->isActive()))
            return;
        
        const std::vector<ot::Block> blkList = pMesh->getLocalBlockList();
        const ot::TreeNode* pNodes = pMesh->getAllElements().data();
        const double current_t= m_uiTimeInfo._m_uiT;
        const unsigned int DOF = m_uiEVar.GetDof();
        
        double current_t_adv=current_t;
        const double dt_finest = m_uiTimeInfo._m_uiTh;
        const double dt_coarset = (1u<<(m_uiLevMax-m_uiLevMin)) * m_uiTimeInfo._m_uiTh;
        
        const unsigned int  finest_t = 1u<<(m_uiLevMax - m_uiLevMax); // finest  time level.
        const unsigned int coarset_t = 1u<<(m_uiLevMax - m_uiLevMin); // coarset time level. 


        assert(blkList.size() > 0);

        for(unsigned int blk=0; blk < blkList.size(); blk++)
        {
            const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
            m_uiBlkTimeLevMinMax[2*bLev +0] = LOOK_UP_TABLE_DEFAULT;
            m_uiBlkTimeLevMinMax[2*bLev +1] = 0;
        }

        for(unsigned int blk=0; blk < blkList.size(); blk++)
        {
            const unsigned int BLK_T = m_uiBVec[blk]._time;
            const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
            
            if(BLK_T  < m_uiBlkTimeLevMinMax[2*bLev +0])
                m_uiBlkTimeLevMinMax[2*bLev + 0] = BLK_T;

            if(BLK_T  > m_uiBlkTimeLevMinMax[2*bLev + 1])
                m_uiBlkTimeLevMinMax[2*bLev + 1] = BLK_T;
        }

        
        for(unsigned int rk=1; rk <= m_uiNumStages; rk++ )
        {
            const unsigned int BLK_S = rk;
            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                const unsigned int BLK_T = m_uiBVec[blk]._time;
                const unsigned int NN   =  m_uiBVec[blk]._vec[0].getSz();
                const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
                const unsigned int BLK_DT = 1u<<(m_uiLevMax - bLev);

                const T blk_dt = m_uiTimeInfo._m_uiTh * BLK_DT;
                
                bool isBlkSynced = true;
                if(m_uiBlkTimeLevMinMax[2*bLev+0] !=  m_uiBlkTimeLevMinMax[2*bLev+1] )
                    isBlkSynced = false;

                // if(!isBlkSynced)                    
                //     std::cout<<"isBlkSynced : "<<isBlkSynced<<" blk min : "<<m_uiBlkTimeLevMinMax[2*bLev+0]<<" blk max "<< m_uiBlkTimeLevMinMax[2*bLev+1]<<std::endl;
                // skip blocks is they are not synced, and the time is equal to the max time. 
                if( (!isBlkSynced && BLK_T == m_uiBlkTimeLevMinMax[2*bLev+1]) || (pt % BLK_DT !=0) )
                    continue;

                //std::cout<<"[NUTS]: pt: "<<pt<<" blk: "<<blk<<" rk : "<<rk<<" step size: "<<BLK_DT<<" time : "<<BLK_T<<std::endl;

                assert(m_uiBVec[blk]._vec[rk-1].isSynced());

                m_uiBVec[blk]._vec[m_uiNumStages+1].vecCopy(m_uiBVec[blk]._vec[0]); // copy the in vector. 
                T* out_ptr = (T*)m_uiBVec[blk]._vec[m_uiNumStages+1].data();
                
                for(unsigned int s=1;  s < BLK_S; s++)
                {
                    const T* stage_ptr = m_uiBVec[blk]._vec[s].data();
                    for(unsigned int n =0; n < DOF*NN; n++){
                        out_ptr[n] += ( blk_dt * m_uiAij[ (BLK_S-1) * m_uiNumStages + (s-1)] * stage_ptr[n]);
                        //printf("s[%d]: %f\n", n, out_ptr[n]);
                        //printf("s: %d : stage[%d]: %f\n", s,  n, stage_ptr[n]);
                    }
                }

                m_uiBVec[blk]._vec[BLK_S].computeVec(m_uiAppCtx, m_uiBVec[blk]._vec[m_uiNumStages+1], current_t + dt_finest*BLK_T);
                m_uiBVec[blk]._vec[BLK_S].mark_unsynced();
                // std::cout<<"rk:"<<rk<<std::endl;
                // m_uiBVec[blk]._vec[BLK_S].dump_vec(std::cout);

                #ifdef __PROFILE_ENUTS__
                    m_uiPt[ENUTSPROFILE::ENUTS_BLK_ZIP].start();
                #endif

                m_uiBVec[blk]._vec[BLK_S].zipDG(pMesh,m_uiStVec[BLK_S-1].GetVecArray(),DOF);
                m_uiBVec[blk]._rks = rk;

                #ifdef __PROFILE_ENUTS__
                    m_uiPt[ENUTSPROFILE::ENUTS_BLK_ZIP].stop();
                #endif

            }

            
            // do the DG vec ghost sync. 
            pMesh->readFromGhostBeginEleDGVec(m_uiStVec[BLK_S-1].GetVecArray(),DOF);
            pMesh->readFromGhostEndEleDGVec(m_uiStVec[BLK_S-1].GetVecArray(),DOF);

            
            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                const unsigned int BLK_T = m_uiBVec[blk]._time;
                const unsigned int NN   =  m_uiBVec[blk]._vec[0].getSz();
                const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
                const unsigned int BLK_DT = 1u<<(m_uiLevMax - bLev);

                bool isBlkSynced = true;
                if(m_uiBlkTimeLevMinMax[2*bLev+0] !=  m_uiBlkTimeLevMinMax[2*bLev+1] )
                    isBlkSynced = false;

                // skip blocks is they are not synced, and the time is equal to the max time. 
                if( (!isBlkSynced && BLK_T == m_uiBlkTimeLevMinMax[2*bLev+1]) || (pt % BLK_DT !=0) )
                    continue;

                sync_blk_timestep(blk,BLK_S);
                m_uiBVec[blk]._vec[BLK_S].mark_synced();
                //std::cout<<"[NUTS]: pt :  "<<pt<<" blk : "<<blk<<" rk : "<<rk<<" step size: "<<BLK_DT<<" is synced: "<<m_uiBVec[blk]._vec[BLK_S].isSynced()<<std::endl;

            }
                
            
        }

        // compute the time step vector and increment time. 
        for(unsigned int blk =0; blk < blkList.size(); blk++)
        {

            const unsigned int BLK_T = m_uiBVec[blk]._time;
            const unsigned int NN   =  m_uiBVec[blk]._vec[0].getSz();
            const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
            const unsigned int BLK_DT = 1u<<(m_uiLevMax - bLev);
            const T blk_dt = m_uiTimeInfo._m_uiTh * BLK_DT;
            
            bool isBlkSynced = true;
            if(m_uiBlkTimeLevMinMax[2*bLev+0] !=  m_uiBlkTimeLevMinMax[2*bLev+1] )
                isBlkSynced = false;

            // skip blocks is they are not synced, and the time is equal to the max time. 
            if( (!isBlkSynced && BLK_T == m_uiBlkTimeLevMinMax[2*bLev+1]) || (pt % BLK_DT !=0) )
                continue;

            T* out_ptr = (T*)m_uiBVec[blk]._vec[0].data();
            for(unsigned int s=1;  s <= m_uiNumStages; s++)
            {
                assert(m_uiBVec[blk]._vec[s].isSynced());
                const T* stage_ptr = m_uiBVec[blk]._vec[s].data();
                for(unsigned int n =0; n < DOF*NN; n++)
                    out_ptr[n] += (blk_dt * m_uiBi[(s-1)]*stage_ptr[n]);
            }

            m_uiBVec[blk]._time += 1u<<(m_uiLevMax -bLev); 
            m_uiBVec[blk]._rks=0;
            m_uiBVec[blk]._vec[0].mark_synced();

            for(unsigned int s=1;  s <= m_uiNumStages; s++)
                m_uiBVec[blk]._vec[s].mark_unsynced();

        }


        update_ele_timestep();
        
        pMesh->readFromGhostBeginElementVec(m_uiEleTime.data(),1);
        pMesh->readFromGhostEndElementVec(m_uiEleTime.data(),1);

        pMesh->waitActive();
        


        
        
    }
   
    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::evolve()
    {

        #ifdef __PROFILE_ENUTS__
            m_uiPt[ENUTSPROFILE::ENUTS_EVOLVE].start();
        #endif

        ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        const double current_t= m_uiTimeInfo._m_uiT;
        double current_t_adv=current_t;
        const double dt_finest = m_uiTimeInfo._m_uiTh;
        const double dt_coarset = (1u<<(m_uiLevMax-m_uiLevMin)) * m_uiTimeInfo._m_uiTh;
        // Assumption: m_uiEvar is time synced accross all the blocks. 

        const unsigned int  finest_t = 1u<<(m_uiLevMax - m_uiLevMax); // finest  time level.
        const unsigned int coarset_t = 1u<<(m_uiLevMax - m_uiLevMin); // coarset time level. 
        

        if(pMesh->isActive())
        {

            const unsigned int rank = pMesh->getMPIRank();
            const unsigned int npes = pMesh->getMPICommSize();
            
            assert (m_uiEVar.GetDof()== m_uiEVecTmp[0].GetDof());
            assert (m_uiEVar.GetDof()== m_uiEVecTmp[1].GetDof());
            
            const unsigned int DOF = m_uiEVar.GetDof();
            MPI_Comm comm = pMesh->getMPICommunicator();

            const ot::TreeNode* const pNodes = pMesh->getAllElements().data();
            m_uiAppCtx->pre_timestep(m_uiEVar);

            #ifdef __PROFILE_ENUTS__
                m_uiPt[ENUTSPROFILE::ENUTS_BLK_UNZIP].start();
            #endif

            m_uiAppCtx->unzip(m_uiEVar,m_uiEvarUzip ,m_uiAppCtx->get_async_batch_sz());  // unzip m_uiEVec to m_uiEVarUnzip
            const std::vector<ot::Block> blkList = pMesh->getLocalBlockList();

            // initialize the block async vectors
            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                m_uiBVec[blk]._time = 0;
                m_uiBVec[blk]._rks  = 0;
                m_uiBVec[blk]._vec[0].copyFromUnzip(pMesh,m_uiEvarUzip.GetVecArray(), true, DOF);
                m_uiBVec[blk]._vec[0].mark_synced();

                for(unsigned int s=1;  s <= m_uiNumStages; s++)
                    m_uiBVec[blk]._vec[s].mark_unsynced();
            }

            #ifdef __PROFILE_ENUTS__
                m_uiPt[ENUTSPROFILE::ENUTS_BLK_UNZIP].stop();
            #endif

            
            for(unsigned int ele = pMesh->getElementPreGhostBegin(); ele < pMesh->getElementPostGhostEnd(); ele ++)
                m_uiEleTime[ele]=0;

            this->update_ele_timestep();

            
            pMesh->readFromGhostBeginElementVec(m_uiEleTime.data(),1);
            pMesh->readFromGhostEndElementVec(m_uiEleTime.data(),1);
            
            
            for(unsigned int pt=0; pt< coarset_t; pt ++)
            {
                
                // unsigned int numBlkEvolved = 0; 
                // for(unsigned int blk =0; blk < blkList.size(); blk++)
                // {
                //     const unsigned int BLK_T = m_uiBVec[blk]._time;
                //     const unsigned int NN   =  m_uiBVec[blk]._vec[0].getSz();
                //     const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
                //     const unsigned int BLK_DT = 1u<<(m_uiLevMax - bLev);
                //     if( pt% BLK_DT !=0 )
                //         continue;
                //     numBlkEvolved++;
                // }

                // std::cout<<"[ENUTS] : pt: "<<pt<<" number of blocks evolved : "<<numBlkEvolved<<" out of "<<blkList.size()<<std::endl;
                this->partial_evolve(pt);
                
            }


            #ifdef __PROFILE_ENUTS__
                m_uiPt[ENUTSPROFILE::ENUTS_BLK_ZIP].start();
            #endif

            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                const unsigned int BLK_T = m_uiBVec[blk]._time;
                const unsigned int NN   =  m_uiBVec[blk]._vec[0].getSz();
                const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
                const unsigned int BLK_DT = 1u<<(m_uiLevMax - bLev);

                m_uiBVec[blk]._vec[0].zip(pMesh, m_uiEVar.GetVecArray(),DOF);

            }

            #ifdef __PROFILE_ENUTS__
                m_uiPt[ENUTSPROFILE::ENUTS_BLK_ZIP].stop();
            #endif
            

        }

        

        m_uiAppCtx->increment_ts_info((T)coarset_t,1);
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        pMesh->waitAll();

        #ifdef __PROFILE_ENUTS__
            m_uiPt[ENUTSPROFILE::ENUTS_EVOLVE].stop();
        #endif


    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::remesh(unsigned int grain_sz, double ld_tol, unsigned int sf_k)
    {
        ot::Mesh* pMesh   = m_uiAppCtx->get_mesh();
        ot::Mesh* newMesh = m_uiAppCtx->remesh(grain_sz,ld_tol,sf_k);
        
        DendroIntL oldElements = pMesh->getNumLocalMeshElements();
        DendroIntL newElements = newMesh->getNumLocalMeshElements();

        DendroIntL oldElements_g, newElements_g;

        par::Mpi_Reduce(&oldElements,&oldElements_g,1,MPI_SUM,0,pMesh->getMPIGlobalCommunicator());
        par::Mpi_Reduce(&newElements,&newElements_g,1,MPI_SUM,0,newMesh->getMPIGlobalCommunicator());

        // if(!(m_uiMesh->getMPIRankGlobal()))
        //     std::cout<<"[LTS]: remesh : "<<m_uiTinfo._m_uiStep<<"\ttime : "<<m_uiTinfo._m_uiT<<"\told mesh: "<<oldElements_g<<"\tnew mesh:"<<newElements_g<<std::endl;

        const unsigned int DOF = m_uiEvarDG.GetDof();
        DVec evarDG;
        evarDG.VecCreateDG(pMesh,true,DOF);

        std::vector<unsigned int> eleTimeVec;
        eleTimeVec.resize(newMesh->getAllElements().size(),0);

        blk_vec_to_zipDG(m_uiBVec.data(),m_uiEvarDG,0);

        pMesh->interGridTransfer_DG(m_uiEvarDG.GetVecArray(), evarDG.GetVecArray(), newMesh, m_uiEvarDG.GetDof());
        pMesh->interGridTransferCellVec(m_uiEleTime.data(),eleTimeVec.data(),newMesh,1,ot::INTERGRID_TRANSFER_MODE::CELLVEC_CPY);
        
        m_uiAppCtx->grid_transfer(newMesh,true,false,false);
        m_uiAppCtx->update_app_vars();

        std::swap(pMesh,newMesh);
        delete newMesh;

        

        m_uiAppCtx->set_mesh(pMesh);
        
        this->deallocate_internal_vars();
        this->free_data_structures();
        
        this->init_data_structures();
        this->allocate_internal_vars();

        // note: I don't think we need a ghost exchange since we are dealing with the DG vectors. 
        // pMesh->readFromGhostBeginEleDGVec(evarDG.GetVecArray(),evarDG.GetDof());
        // pMesh->readFromGhostEndEleDGVec(evarDG.GetVecArray(),evarDG.GetDof());

        vecDG_to_blk_vec(evarDG,m_uiBVec.data(),0);        

        std::swap(m_uiEleTime,eleTimeVec);
        eleTimeVec.clear();
        

        return;
        

    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::evolve_with_remesh(double sync_time, unsigned int remesh_freq)
    {

        #ifdef __PROFILE_ENUTS__
            m_uiPt[ENUTSPROFILE::ENUTS_EVOLVE].start();
        #endif

        ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        const double current_t= m_uiTimeInfo._m_uiT;
        double current_t_adv=current_t;
        
        const double dt_finest = m_uiTimeInfo._m_uiTh;
        const double dt_coarset = (1u<<(m_uiLevMax-m_uiLevMin)) * m_uiTimeInfo._m_uiTh;

        // Assumption: m_uiEvar is time synced accross all the blocks. 
        const unsigned int  finest_t = 1u<<(m_uiLevMax - m_uiLevMax); // finest  time level.
        const unsigned int coarset_t = 1u<<(m_uiLevMax - m_uiLevMin); // coarset time level. 
        const unsigned int sync_level = m_uiLevMin;

        


        if(pMesh->isActive())
        {

            const unsigned int rank = pMesh->getMPIRank();
            const unsigned int npes = pMesh->getMPICommSize();
            
            assert (m_uiEVar.GetDof()== m_uiEVecTmp[0].GetDof());
            assert (m_uiEVar.GetDof()== m_uiEVecTmp[1].GetDof());
            
            const unsigned int DOF = m_uiEVar.GetDof();
            MPI_Comm comm = pMesh->getMPICommunicator();

            const ot::TreeNode* const pNodes = pMesh->getAllElements().data();
            const std::vector<ot::Block> blkList = pMesh->getLocalBlockList();
            
            m_uiAppCtx->pre_timestep(m_uiEVar);

            #ifdef __PROFILE_ENUTS__
                m_uiPt[ENUTSPROFILE::ENUTS_BLK_UNZIP].start();
            #endif

            m_uiAppCtx->unzip(m_uiEVar,m_uiEvarUzip ,m_uiAppCtx->get_async_batch_sz());  // unzip m_uiEVec to m_uiEVarUnzip
            

            // initialize the block async vectors
            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                m_uiBVec[blk]._time = 0;
                m_uiBVec[blk]._rks  = 0;
                m_uiBVec[blk]._vec[0].copyFromUnzip(pMesh,m_uiEvarUzip.GetVecArray(), true, DOF);
                m_uiBVec[blk]._vec[0].mark_synced();

                for(unsigned int s=1;  s <= m_uiNumStages; s++)
                    m_uiBVec[blk]._vec[s].mark_unsynced();
            }

            #ifdef __PROFILE_ENUTS__
                m_uiPt[ENUTSPROFILE::ENUTS_BLK_UNZIP].stop();
            #endif

            
            for(unsigned int ele = pMesh->getElementPreGhostBegin(); ele < pMesh->getElementPostGhostEnd(); ele ++)
                m_uiEleTime[ele]=0;

            this->update_ele_timestep();

            
            pMesh->readFromGhostBeginElementVec(m_uiEleTime.data(),1);
            pMesh->readFromGhostEndElementVec(m_uiEleTime.data(),1);


        }


        unsigned int pt=0;
        std::vector<unsigned int> blk_time_min_max;
        blk_time_min_max.resize(m_uiMaxDepth,0);

        while( pt < coarset_t)
        {

            if( (pt % remesh_freq) == 0)
            {
                bool isRemesh =false ; // pMesh->isReMeshUnzip()

                if(isRemesh)
                {
                    // create new mesh. 
                    // do intergrid trasfer 
                        // 1. zipDG vector. 
                        // 2. timeCell vector
                    
                    // compute new min and max levels. 
                    this->remesh();
                }

            }

            pMesh = m_uiAppCtx->get_mesh();
            if(pMesh->isActive())
                this->partial_evolve(pt);

            pt++;
            

        }

        pMesh = m_uiAppCtx->get_mesh();

        if(pMesh->isActive())
        {   
            
            const std::vector<ot::Block> blkList = pMesh->getLocalBlockList();

            #ifdef __PROFILE_ENUTS__
                m_uiPt[ENUTSPROFILE::ENUTS_BLK_ZIP].start();
            #endif

            const unsigned int DOF = m_uiEVar.GetDof();
            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                //const unsigned int BLK_T = m_uiBVec[blk]._time;
                //const unsigned int NN   =  m_uiBVec[blk]._vec[0].getSz();
                //const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
                //const unsigned int BLK_DT = 1u<<(m_uiLevMax - bLev);
                m_uiBVec[blk]._vec[0].zip(pMesh, m_uiEVar.GetVecArray(),DOF);
            }

            #ifdef __PROFILE_ENUTS__
                m_uiPt[ENUTSPROFILE::ENUTS_BLK_ZIP].stop();
            #endif



        }

        
        
        
        m_uiAppCtx->increment_ts_info((T)coarset_t,1);
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        pMesh->waitAll();

        #ifdef __PROFILE_ENUTS__
            m_uiPt[ENUTSPROFILE::ENUTS_EVOLVE].stop();
        #endif


    }


}


