//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Header file for the GR simulation.
*/
//

#include "quadgrav.h"
#include "quadgravUtils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "rk4quadgrav.h"
#include "octUtils.h"

int main (int argc, char** argv)
{


    if(argc<2)
        std::cout<<"Usage: "<<argv[0]<<" paramFile"<<std::endl;

    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    quadgrav::timer::initFlops();

    quadgrav::timer::total_runtime.start();


    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    quadgrav::readParamFile(argv[1],comm);



    if(rank==1|| npes==1)
    {
        std::cout<<"parameters read: "<<std::endl;

        std::cout<<YLW<<"\tnpes :"<<npes<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_DIM :"<<quadgrav::QUADGRAV_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_IO_OUTPUT_FREQ :"<<quadgrav::QUADGRAV_IO_OUTPUT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_REMESH_TEST_FREQ :"<<quadgrav::QUADGRAV_REMESH_TEST_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_CHECKPT_FREQ :"<<quadgrav::QUADGRAV_CHECKPT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_RESTORE_SOLVER :"<<quadgrav::QUADGRAV_RESTORE_SOLVER<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ENABLE_BLOCK_ADAPTIVITY :"<<quadgrav::QUADGRAV_ENABLE_BLOCK_ADAPTIVITY<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_VTU_FILE_PREFIX :"<<quadgrav::QUADGRAV_VTU_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_CHKPT_FILE_PREFIX :"<<quadgrav::QUADGRAV_CHKPT_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_PROFILE_FILE_PREFIX :"<<quadgrav::QUADGRAV_PROFILE_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_IO_OUTPUT_GAP :"<<quadgrav::QUADGRAV_IO_OUTPUT_GAP<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_DENDRO_GRAIN_SZ :"<<quadgrav::QUADGRAV_DENDRO_GRAIN_SZ<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ASYNC_COMM_K :"<<quadgrav::QUADGRAV_ASYNC_COMM_K<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_DENDRO_AMR_FAC :"<<quadgrav::QUADGRAV_DENDRO_AMR_FAC<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_CFL_FACTOR:"<<quadgrav::QUADGRAV_CFL_FACTOR<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_WAVELET_TOL :"<<quadgrav::QUADGRAV_WAVELET_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_LOAD_IMB_TOL :"<<quadgrav::QUADGRAV_LOAD_IMB_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_RK45_TIME_BEGIN :"<<quadgrav::QUADGRAV_RK45_TIME_BEGIN<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_RK45_TIME_END :"<<quadgrav::QUADGRAV_RK45_TIME_END<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_RK45_TIME_STEP_SIZE :"<<quadgrav::QUADGRAV_RK45_TIME_STEP_SIZE<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_RK45_DESIRED_TOL :"<<quadgrav::QUADGRAV_RK45_DESIRED_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_COMPD_MIN : ( :"<<quadgrav::QUADGRAV_COMPD_MIN[0]<<" ,"<<quadgrav::QUADGRAV_COMPD_MIN[1]<<","<<quadgrav::QUADGRAV_COMPD_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_COMPD_MAX : ( :"<<quadgrav::QUADGRAV_COMPD_MAX[0]<<" ,"<<quadgrav::QUADGRAV_COMPD_MAX[1]<<","<<quadgrav::QUADGRAV_COMPD_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_BLK_MIN : ( :"<<quadgrav::QUADGRAV_BLK_MIN_X<<" ,"<<quadgrav::QUADGRAV_BLK_MIN_Y<<","<<quadgrav::QUADGRAV_BLK_MIN_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_BLK_MAX : ( :"<<quadgrav::QUADGRAV_BLK_MAX_X<<" ,"<<quadgrav::QUADGRAV_BLK_MAX_Y<<","<<quadgrav::QUADGRAV_BLK_MAX_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_OCTREE_MIN : ( :"<<quadgrav::QUADGRAV_OCTREE_MIN[0]<<" ,"<<quadgrav::QUADGRAV_OCTREE_MIN[1]<<","<<quadgrav::QUADGRAV_OCTREE_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_OCTREE_MAX : ( :"<<quadgrav::QUADGRAV_OCTREE_MAX[0]<<" ,"<<quadgrav::QUADGRAV_OCTREE_MAX[1]<<","<<quadgrav::QUADGRAV_OCTREE_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tKO_DISS_SIGMA :"<<quadgrav::KO_DISS_SIGMA<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_TYPE:"<<quadgrav::QUADGRAV_ID_TYPE<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_AMP1:"<<quadgrav::QUADGRAV_ID_AMP1<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_AMP2:"<<quadgrav::QUADGRAV_ID_AMP2<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_DELTA1:"<<quadgrav::QUADGRAV_ID_DELTA1<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_DELTA2:"<<quadgrav::QUADGRAV_ID_DELTA2<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_XC1:"<<quadgrav::QUADGRAV_ID_XC1<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_YC1:"<<quadgrav::QUADGRAV_ID_YC1<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_ZC1:"<<quadgrav::QUADGRAV_ID_ZC1<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_XC2:"<<quadgrav::QUADGRAV_ID_XC2<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_YC2:"<<quadgrav::QUADGRAV_ID_YC2<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_ZC2:"<<quadgrav::QUADGRAV_ID_ZC2<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_EPSX1:"<<quadgrav::QUADGRAV_ID_EPSX1<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_EPSY1:"<<quadgrav::QUADGRAV_ID_EPSY1<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_EPSZ1:"<<quadgrav::QUADGRAV_ID_EPSY1<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_EPSX2:"<<quadgrav::QUADGRAV_ID_EPSX2<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_EPSY2:"<<quadgrav::QUADGRAV_ID_EPSY2<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_EPSZ2:"<<quadgrav::QUADGRAV_ID_EPSY2<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_R1:"<<quadgrav::QUADGRAV_ID_R1<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_R2:"<<quadgrav::QUADGRAV_ID_R2<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_NU1:"<<quadgrav::QUADGRAV_ID_NU1<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_NU2:"<<quadgrav::QUADGRAV_ID_NU2<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ID_OMEGA:"<<quadgrav::QUADGRAV_ID_OMEGA<<NRM<<std::endl;
        
        

        std::cout<<YLW<<"\tQUADGRAV_DIM :"<<quadgrav::QUADGRAV_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_MAXDEPTH :"<<quadgrav::QUADGRAV_MAXDEPTH<<NRM<<std::endl;

        std::cout<<YLW<<"\tQUADGRAV_NUM_REFINE_VARS :"<<quadgrav::QUADGRAV_NUM_REFINE_VARS<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_REFINE_VARIABLE_INDICES :[";
        for(unsigned int i=0;i<quadgrav::QUADGRAV_NUM_REFINE_VARS-1;i++)
            std::cout<<quadgrav::QUADGRAV_REFINE_VARIABLE_INDICES[i]<<", ";
        std::cout<<quadgrav::QUADGRAV_REFINE_VARIABLE_INDICES[quadgrav::QUADGRAV_NUM_REFINE_VARS-1]<<"]"<<NRM<<std::endl;

        std::cout<<YLW<<"\tQUADGRAV_NUM_EVOL_VARS_VTU_OUTPUT :"<<quadgrav::QUADGRAV_NUM_EVOL_VARS_VTU_OUTPUT<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_VTU_OUTPUT_EVOL_INDICES :[";
        for(unsigned int i=0;i<quadgrav::QUADGRAV_NUM_EVOL_VARS_VTU_OUTPUT-1;i++)
            std::cout<<quadgrav::QUADGRAV_VTU_OUTPUT_EVOL_INDICES[i]<<", ";
        std::cout<<quadgrav::QUADGRAV_VTU_OUTPUT_EVOL_INDICES[quadgrav::QUADGRAV_NUM_EVOL_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;



    }

    _InitializeHcurve(quadgrav::QUADGRAV_DIM);
    m_uiMaxDepth=quadgrav::QUADGRAV_MAXDEPTH;

    if(quadgrav::QUADGRAV_NUM_VARS%quadgrav::QUADGRAV_ASYNC_COMM_K!=0)
    {
        if(!rank) std::cout<<"[overlap communication error]: total QUADGRAV_NUM_VARS: "<<quadgrav::QUADGRAV_NUM_VARS<<" is not divisable by QUADGRAV_ASYNC_COMM_K: "<<quadgrav::QUADGRAV_ASYNC_COMM_K<<std::endl;
        exit(0);
    }

    quadgrav::QUADGRAV_RK45_TIME_STEP_SIZE=quadgrav::QUADGRAV_CFL_FACTOR*(quadgrav::QUADGRAV_COMPD_MAX[0]-quadgrav::QUADGRAV_COMPD_MIN[0])*(1.0/(double)(1u<<quadgrav::QUADGRAV_MAXDEPTH));

    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){quadgrav::initData(x,y,z,var);};
    std::function<void(double,double,double,double,double*)> u_x_t=[](double x,double y,double z,double t,double*var){quadgrav::analyticalSol(x,y,z,t,var);};
    //std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){quadgrav::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars=quadgrav::QUADGRAV_NUM_VARS;
    unsigned int varIndex[interpVars];
    for(unsigned int i=0;i<quadgrav::QUADGRAV_NUM_VARS;i++)
        varIndex[i]=i;

    
    DendroIntL localSz,globalSz;
    double t_stat;
    double t_stat_g[3];

    quadgrav::timer::t_f2o.start();

    if(quadgrav::QUADGRAV_ENABLE_BLOCK_ADAPTIVITY)
    {
        if(!rank) std::cout<<YLW<<"Using block adaptive mesh. AMR disabled "<<NRM<<std::endl;
        const Point pt_min(quadgrav::QUADGRAV_BLK_MIN_X,quadgrav::QUADGRAV_BLK_MIN_Y,quadgrav::QUADGRAV_BLK_MIN_Z);
        const Point pt_max(quadgrav::QUADGRAV_BLK_MAX_X,quadgrav::QUADGRAV_BLK_MAX_Y,quadgrav::QUADGRAV_BLK_MAX_Z);

        quadgrav::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-2,m_uiMaxDepth,comm);
    }else
    {

        if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
        function2Octree(f_init,quadgrav::QUADGRAV_NUM_VARS,quadgrav::QUADGRAV_REFINE_VARIABLE_INDICES,quadgrav::QUADGRAV_NUM_REFINE_VARS,tmpNodes,m_uiMaxDepth,quadgrav::QUADGRAV_WAVELET_TOL,quadgrav::QUADGRAV_ELE_ORDER,comm);
        //std::cout<<"f2o else end"<<std::endl;

    }

    quadgrav::timer::t_f2o.stop();

    t_stat=quadgrav::timer::t_f2o.seconds;
    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    localSz=tmpNodes.size();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);

    if(!rank) std::cout<<GRN<<" function to octree max (s): "<<t_stat_g[2]<<NRM<<std::endl;
    if(!rank) std::cout<<GRN<<" function to octree # octants : "<<globalSz<<NRM<<std::endl;

    par::Mpi_Bcast(&globalSz,1,0,comm);
    const unsigned int grainSz=quadgrav::QUADGRAV_DENDRO_GRAIN_SZ;//DENDRO_DEFAULT_GRAIN_SZ;

    bool isActive;
    MPI_Comm commActive;
    const int p_npes_prev=binOp::getPrevHighestPowerOfTwo((globalSz/grainSz));
    const int p_npes_next=binOp::getNextHighestPowerOfTwo((globalSz/grainSz));

    int p_npes=globalSz/grainSz;
    (std::abs(p_npes_prev-p_npes)<=std::abs(p_npes_next-p_npes)) ? p_npes=p_npes_prev : p_npes=p_npes_next;

    if(p_npes>npes) p_npes=npes;
    // quick fix to enforce the npes>=2 for any given grain size.
    if(p_npes<=1 && npes>1) p_npes=2;

    if(p_npes==npes)
    {
        MPI_Comm_dup(comm,&commActive);
        isActive=true;

    }else
    {
        //isActive=(rank*grainSz<globalSz);
        isActive=isRankSelected(npes,rank,p_npes);
        par::splitComm2way(isActive,&commActive,comm);

    }

    shrinkOrExpandOctree(tmpNodes,quadgrav::QUADGRAV_LOAD_IMB_TOL,DENDRO_DEFAULT_SF_K,isActive,commActive,comm);

    if(!isActive)
        if(tmpNodes.size()!=0)
            std::cout<<" rank_g: "<<rank<<" isActive: "<<isActive<<" f2O octants: "<<tmpNodes.size()<<std::endl;




    std::vector<ot::TreeNode> balOct;
    localSz=0;
    if(isActive)
    {

        int rank_active,npes_active;

        MPI_Comm_size(commActive,&npes_active);
        MPI_Comm_rank(commActive,&rank_active);

        if(!rank_active) std::cout<<"[MPI_COMM_SWITCH]: "<<npes_active<<std::endl;

        ot::TreeNode root(quadgrav::QUADGRAV_DIM,quadgrav::QUADGRAV_MAXDEPTH);
        std::vector<ot::TreeNode> tmpVec;
        quadgrav::timer::t_cons.start();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,quadgrav::QUADGRAV_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,quadgrav::QUADGRAV_SPLIT_FIX,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,quadgrav::QUADGRAV_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,quadgrav::QUADGRAV_SPLIT_FIX,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        quadgrav::timer::t_cons.stop();
        t_stat=quadgrav::timer::t_cons.seconds;

        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
        t_stat_g[1]=t_stat_g[1]/(double)rank_active;

        localSz=tmpNodes.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,commActive);

        if(!rank_active) std::cout<<GRN<<"remove duplicates + octree construction (s): "<<t_stat_g[2]<<NRM<<std::endl;
        if(!rank_active) std::cout<<GRN<<" # const. octants: "<<globalSz<<NRM<<std::endl;


        quadgrav::timer::t_bal.start();

        SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,quadgrav::QUADGRAV_LOAD_IMB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,quadgrav::QUADGRAV_SPLIT_FIX,commActive);
        tmpNodes.clear();

        quadgrav::timer::t_bal.stop();

        t_stat=quadgrav::timer::t_bal.seconds;
        par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,commActive);
        par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,commActive);
        t_stat_g[1]=t_stat_g[1]/(double)rank_active;

        if(!rank_active) std::cout<<GRN<<" 2:1 balancing max (s): "<<t_stat_g[2]<<NRM<<std::endl;
        localSz=balOct.size();


    }
    MPI_Comm_free(&commActive);

    // all reduce act as barrier to sync all procs.
    par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
    if(!rank) std::cout<<GRN<<" balanced # octants : "<<globalSz<<NRM<<std::endl;

    quadgrav::timer::t_mesh.start();

    ot::Mesh * mesh=new ot::Mesh(balOct,1,quadgrav::QUADGRAV_ELE_ORDER,comm,true,ot::SM_TYPE::FDM,quadgrav::QUADGRAV_DENDRO_GRAIN_SZ,quadgrav::QUADGRAV_LOAD_IMB_TOL,quadgrav::QUADGRAV_SPLIT_FIX);

    quadgrav::timer::t_mesh.stop();

    t_stat=quadgrav::timer::t_mesh.seconds;
    par::Mpi_Reduce(&t_stat,t_stat_g,1,MPI_MIN,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+1,1,MPI_SUM,0,comm);
    par::Mpi_Reduce(&t_stat,t_stat_g+2,1,MPI_MAX,0,comm);
    t_stat_g[1]=t_stat_g[1]/(double)npes;

    localSz=mesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
    if(!rank) std::cout<<GRN<<" # of CG nodes (vertices) : "<<globalSz<<NRM<<std::endl;
    if(!rank)
    {
        std::cout<< GRN<<"Mesh generation time (max): "<<t_stat_g[2]<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" e2e (min,mean,max): "<<"( "<<t_e2e_g[0]<<"\t"<<t_e2e_g[1]<<"\t"<<t_e2e_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" e2n (min,mean,max): "<<"( "<<t_e2n_g[0]<<"\t"<<t_e2n_g[1]<<"\t"<<t_e2n_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" sm (min,mean,max): "<<"( "<<t_sm_g[0]<<"\t"<<t_sm_g[1]<<"\t"<<t_sm_g[2]<<" )"<<NRM<<std::endl;
        std::cout<<"\t"<<GRN<<" blk (min,mean,max): "<<"( "<<t_blk_g[0]<<"\t"<<t_blk_g[1]<<"\t"<<t_blk_g[2]<<" )"<<NRM<<std::endl;
    }


    ode::solver::RK4_QUADGRAV rk_quadgrav(mesh,quadgrav::QUADGRAV_RK45_TIME_BEGIN,quadgrav::QUADGRAV_RK45_TIME_END,quadgrav::QUADGRAV_RK45_TIME_STEP_SIZE);
    //ode::solver::RK3_QUADGRAV rk_quadgrav(mesh,quadgrav::QUADGRAV_RK45_TIME_BEGIN,quadgrav::QUADGRAV_RK45_TIME_END,quadgrav::QUADGRAV_RK45_TIME_STEP_SIZE);

    if(quadgrav::QUADGRAV_RESTORE_SOLVER==1)
        rk_quadgrav.restoreCheckPoint(quadgrav::QUADGRAV_CHKPT_FILE_PREFIX.c_str(),comm);

    quadgrav::timer::t_rkSolve.start();
    rk_quadgrav.rkSolve();
    quadgrav::timer::t_rkSolve.stop();

    quadgrav::timer::total_runtime.stop();
    rk_quadgrav.freeMesh();
    //quadgrav::timer::profileInfo(quadgrav::QUADGRAV_PROFILE_FILE_PREFIX.c_str(),mesh);
    //delete mesh;
    MPI_Finalize();

    return 0;
}
