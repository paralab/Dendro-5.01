//
// Created by milinda on 3/31/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains code to apply a give stencil for a given octree.
*/
//

#include "TreeNode.h"
#include "mpi.h"
#include "genPts_par.h"
#include "sfcSort.h"
#include "mesh.h"
#include "dendro.h"
#include "dendroIO.h"
#include "octUtils.h"
#include "functional"
#include "fdCoefficient.h"
#include "stencil.h"
#include "oct2vtk.h"
#include "rkTransportUtils.h"
#include "meshTestUtils.h"
#include "rawIO.h"



int main (int argc, char** argv) {


    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if(argc<4)
    {
        if(!rank) std::cout<<"Usage: "<<argv[0]<<" maxDepth wavelet_tol partition_tol eleOrder"<<std::endl;
        return 0;
    }

    m_uiMaxDepth=atoi(argv[1]);
    double wavelet_tol=atof(argv[2]);
    double partition_tol=atof(argv[3]);
    unsigned int eOrder=atoi(argv[4]);

    if(!rank)
    {
        std::cout<<YLW<<"maxDepth: "<<m_uiMaxDepth<<NRM<<std::endl;
        std::cout<<YLW<<"wavelet_tol: "<<wavelet_tol<<NRM<<std::endl;
        std::cout<<YLW<<"partition_tol: "<<partition_tol<<NRM<<std::endl;
        std::cout<<YLW<<"eleOrder: "<<eOrder<<NRM<<std::endl;

    }

    _InitializeHcurve(m_uiDim);

    // function that we need to interpolate.
    double d_min=-0.5;
    double d_max=0.5;
    double dMin[]={-0.5,-0.5,-0.5};
    double dMax[]={0.5,0.5,0.5};

    Point pt_min(-0.5,-0.5,-0.5);
    Point pt_max(0.5,0.5,0.5);
    //@note that based on how the functions are defined (f(x), dxf(x), etc) the compuatational domain is equivalent to the grid domain.
    std::function<double(double,double,double)> func =[d_min,d_max](const double x,const double y,const double z){ return (sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
    std::function<double(double,double,double)> dx_func=[d_min,d_max](const double x,const double y,const double z){ return (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(cos(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
    std::function<double(double,double,double)> dy_func=[d_min,d_max](const double x,const double y,const double z){ return (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*cos(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};
    std::function<double(double,double,double)> dz_func=[d_min,d_max](const double x,const double y,const double z){ return (2*M_PI*(1.0/(1u<<m_uiMaxDepth)*(d_max-d_min)))*(sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*cos(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};


    //std::function<double(double,double,double)> func =[d_min,d_max](const double x,const double y,const double z){ return (sin(2*M_PI*((x/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((y/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min))*sin(2*M_PI*((z/(1u<<m_uiMaxDepth))*(d_max-d_min)+d_min)));};

    /*std::function<double(double,double,double)> func =[](const double x,const double  y,const double z){ return (x*x+y*y+z*z);};
    std::function<double(double,double,double)> dx_func =[](const double x,const double  y,const double z){ return (2*x);};
    std::function<double(double,double,double)> dy_func =[](const double x,const double  y,const double z){ return (2*y);};
    std::function<double(double,double,double)> dz_func =[](const double x,const double  y,const double z){ return (2*z);};*/
    //std::function<double(double,double,double)> func =[](const double x,const double  y,const double z){ return (x+y+z);};
    std::vector<ot::TreeNode> tmpNodes;

    function2Octree(func,tmpNodes,m_uiMaxDepth,wavelet_tol,eOrder,comm);

    DendroIntL localSz,globalSz;
    localSz=tmpNodes.size();
    par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
    if(!rank) std::cout<<GRN<<" # f20. octants: "<<globalSz<<NRM<<std::endl;

    par::Mpi_Bcast(&globalSz,1,0,comm);
    const unsigned int grainSz=100;//DENDRO_DEFAULT_GRAIN_SZ;

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

    shrinkOrExpandOctree(tmpNodes,DENDRO_DEFAULT_LB_TOL,DENDRO_DEFAULT_SF_K,isActive,commActive,comm);

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

        ot::TreeNode root(m_uiDim,m_uiMaxDepth);
        std::vector<ot::TreeNode> tmpVec;


        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,DENDRO_DEFAULT_LB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_REMOVE_DUPLICATES,DENDRO_DEFAULT_SF_K,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();

        SFC::parSort::SFC_treeSort(tmpNodes,tmpVec,tmpVec,tmpVec,DENDRO_DEFAULT_LB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_CONSTRUCT_OCTREE,DENDRO_DEFAULT_SF_K,commActive);
        std::swap(tmpNodes,tmpVec);
        tmpVec.clear();


        localSz=tmpNodes.size();
        par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,commActive);

        if(!rank_active) std::cout<<GRN<<" # const. octants: "<<globalSz<<NRM<<std::endl;




        SFC::parSort::SFC_treeSort(tmpNodes,balOct,balOct,balOct,DENDRO_DEFAULT_LB_TOL,m_uiMaxDepth,root,ROOT_ROTATION,1,TS_BALANCE_OCTREE,DENDRO_DEFAULT_SF_K,commActive);
        tmpNodes.clear();

        localSz=balOct.size();


    }
    MPI_Comm_free(&commActive);

    // all reduce act as barrier to sync all procs.
    par::Mpi_Allreduce(&localSz,&globalSz,1,MPI_SUM,comm);
    if(!rank) std::cout<<GRN<<" balanced # octants : "<<globalSz<<NRM<<std::endl;

    ot::Mesh* mesh=new ot::Mesh(balOct,1,eOrder,comm);


    DendroIntL numVertices=mesh->getNumLocalMeshNodes();
    DendroIntL numVertices_g;

    par::Mpi_Reduce(&numVertices,&numVertices_g,1,MPI_SUM,0,comm);

    if(!rank)
        std::cout<<"Total number of vertices: "<<numVertices_g<<std::endl;


    MPI_Barrier(MPI_COMM_WORLD);
    if(!rank) std::cout<<"mesh generation ended "<<std::endl;
    //io::vtk::mesh2vtu(&mesh, "balOct",0,NULL,NULL,0,NULL,NULL);

    //if(!rank) std::cout<<"fabs test: "<<fabs(1-2)<<std::endl;
    unsigned int nLocalBegin=mesh->getNodeLocalBegin();
    unsigned int nLocalEnd=mesh->getNodeLocalEnd();
    std::vector<double> funcVal;
    std::vector<double> funcValUnZip;
    std::vector<double> dx_funcVal;
    std::vector<double> dx_funcVal1;

    std::vector<double> dy_funcVal;
    std::vector<double> dy_funcVal1;

    std::vector<double> dz_funcVal;
    std::vector<double> dz_funcVal1;

    mesh->createVector(funcVal,func);
    mesh->createUnZippedVector(funcValUnZip);
    mesh->performGhostExchange(funcVal);

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        Nodal Value Check Begin     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    //assert(ot::test::isElementalNodalValuesValid(&mesh,&(*(funcVal.begin())),func,1e-3));
    ot::test::isElementalNodalValuesValid(mesh,&(*(funcVal.begin())),func,1e-3);
    

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        Nodal Value Check end     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }




    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        Nodal to Hanging Computation Check Begin     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    //assert(ot::test::isElementalNodalValuesValid(&mesh,&(*(funcVal.begin())),func,1e-3));
    ot::test::isElementalContributionValid(mesh,func,func,1e-3);

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        Nodal to Hanging Computation Check End     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }



    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        UnZip Test Check Begin     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    mesh->unzip(&(*(funcVal.begin())),&(*(funcValUnZip.begin())));
    //assert(ot::test::isUnzipValid(&mesh,&(*(funcValUnZip.begin())),func,1e-3));
    ot::test::isUnzipValid(mesh,&(*(funcValUnZip.begin())),func,1e-3);

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        UnZip Test Check End     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        interpolation to sphere check begin    "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }
    
    //ot::test::isInterpToSphereValid(mesh,&(*(funcVal.begin())),func,0.1,dMin,dMax,1e-3);
    ot::test::isSphereInterpValid(mesh,&(*(funcVal.begin())),func,0.1,1e-3,pt_min,pt_max);

    //bool isSphereInterpValid(ot::Mesh* pMesh, T* vec, std::function< double(double,double,double) > func, double r, double tol, Point d_min, Point d_max)

    if(!rank)
    {
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;
        std::cout<<"        interpolation to sphere check end     "<<std::endl;
        std::cout<<YLW<<"================================================================================================"<<NRM<<std::endl;

    }



    mesh->createVector(dx_funcVal1,dx_func);
    mesh->createVector(dy_funcVal1,dy_func);
    mesh->createVector(dz_funcVal1,dz_func);

    mesh->createVector(dx_funcVal);
    mesh->createVector(dy_funcVal);
    mesh->createVector(dz_funcVal);


    Stencil<double,5,2> D1_Order4StencilCentered_x(fd::D1_ORDER_4_CENTERED,5,StencilDirection::STENCIL_DIR_X);
    Stencil<double,5,4> D1_Order4StencilBackward_x(fd::D1_ORDER_4_BACKWARD,5,StencilDirection::STENCIL_DIR_X);
    Stencil<double,5,0> D1_Order4StencilForward_x(fd::D1_ORDER_4_FORWARD,5,StencilDirection::STENCIL_DIR_X);

    Stencil<double,5,1> D1_Order4StencilUpWind_x(fd::D1_ORDER_4_UPWIND,5,StencilDirection::STENCIL_DIR_X);
    Stencil<double,5,3> D1_Order4StencilDownWind_x(fd::D1_ORDER_4_DOWNWIND,5,StencilDirection::STENCIL_DIR_X);


    Stencil<double,5,2> D1_Order4StencilCentered_y(fd::D1_ORDER_4_CENTERED,5,StencilDirection::STENCIL_DIR_Y);
    Stencil<double,5,4> D1_Order4StencilBackward_y(fd::D1_ORDER_4_BACKWARD,5,StencilDirection::STENCIL_DIR_Y);
    Stencil<double,5,0> D1_Order4StencilForward_y(fd::D1_ORDER_4_FORWARD,5,StencilDirection::STENCIL_DIR_Y);

    Stencil<double,5,1> D1_Order4StencilUpWind_y(fd::D1_ORDER_4_UPWIND,5,StencilDirection::STENCIL_DIR_Y);
    Stencil<double,5,3> D1_Order4StencilDownWind_y(fd::D1_ORDER_4_DOWNWIND,5,StencilDirection::STENCIL_DIR_Y);

    Stencil<double,5,2> D1_Order4StencilCentered_z(fd::D1_ORDER_4_CENTERED,5,StencilDirection::STENCIL_DIR_Z);
    Stencil<double,5,4> D1_Order4StencilBackward_z(fd::D1_ORDER_4_BACKWARD,5,StencilDirection::STENCIL_DIR_Z);
    Stencil<double,5,0> D1_Order4StencilForward_z(fd::D1_ORDER_4_FORWARD,5,StencilDirection::STENCIL_DIR_Z);

    Stencil<double,5,1> D1_Order4StencilUpWind_z(fd::D1_ORDER_4_UPWIND,5,StencilDirection::STENCIL_DIR_Z);
    Stencil<double,5,3> D1_Order4StencilDownWind_z(fd::D1_ORDER_4_DOWNWIND,5,StencilDirection::STENCIL_DIR_Z);



    mesh->performGhostExchange(funcVal);



#if 0

    double* vars[3];
    vars[0]=&(*(funcVal.begin()));
    vars[1]=&(*(dy_funcVal1.begin()));
    vars[2]=&(*(dz_funcVal1.begin()));

    std::cout<<"write begin: "<<std::endl;
    io::varToRawData((const ot::Mesh*)mesh,(const double **)vars,3,NULL,"dendro_raw");
    std::cout<<"write end: "<<std::endl;
#endif
   /* std::vector<double> unZipIn;
    std::vector<double> unZipOut;
    //@milinda : If you are using grad functions defined in rkTransport utills, make sure to use the correct domain paramerters when computing the h.
    mesh.createUnZippedVector(unZipIn,0.0);
    mesh.createUnZippedVector(unZipOut,0.0);

    mesh.unzip(&(*(funcVal.begin())),&(*(unZipIn.begin())));

    grad(&mesh,0,&(*(unZipIn.begin())),&(*(unZipOut.begin())));
    mesh.zip(&(*(unZipOut.begin())),&(*(dx_funcVal.begin())));

    grad(&mesh,1,&(*(unZipIn.begin())),&(*(unZipOut.begin())));
    mesh.zip(&(*(unZipOut.begin())),&(*(dy_funcVal.begin())));

    grad(&mesh,2,&(*(unZipIn.begin())),&(*(unZipOut.begin())));
    mesh.zip(&(*(unZipOut.begin())),&(*(dz_funcVal.begin())));

    //std::cout<<"zip vec: "<<funcVal.size()<<" unZipVec: "<<unZipIn.size()<<std::endl;

    unZipOut.clear();
    unZipIn.clear();*/

    const std::vector<ot::Block> blockList=mesh->getLocalBlockList();
    std::vector<ot::TreeNode> localBlocks;

    for(unsigned int e=0;e<blockList.size();e++)
        localBlocks.push_back(blockList[e].getBlockNode());

    treeNodesTovtk(localBlocks,rank,"blocks");

    mesh->applyStencil(funcVal, dx_funcVal, D1_Order4StencilCentered_x, D1_Order4StencilBackward_x, D1_Order4StencilForward_x);
    mesh->applyStencil(funcVal, dy_funcVal, D1_Order4StencilCentered_y, D1_Order4StencilBackward_y, D1_Order4StencilForward_y);
    mesh->applyStencil(funcVal, dz_funcVal, D1_Order4StencilCentered_z, D1_Order4StencilBackward_z, D1_Order4StencilForward_z);

    // Need to do the gost exchange before wirting vtu files.
    mesh->performGhostExchange(funcVal);
    mesh->performGhostExchange(dx_funcVal1);
    mesh->performGhostExchange(dx_funcVal);

    unsigned int numPVars=3;
    const char * pVarNames [] ={"F(x)","dxF(x)","SxF(x)"};
    double * pVarVal []={&(*(funcVal.begin())),&(*(dx_funcVal1.begin())),&(*(dx_funcVal.begin()))};

    const char *fVarNames []={"Time","Cycle"};
    double fVarVal []={0.1,1};

    io::vtk::mesh2vtuFine(mesh, "f_val",2,(const char **)fVarNames,(const double *)fVarVal,numPVars,(const char** )pVarNames,(const double **)pVarVal);



    /*mesh.vectorToVTK(funcVal,"f");
    mesh.vectorToVTK(dx_funcVal1,"dx_f");
    mesh.vectorToVTK(dx_funcVal,"dx_f_approx");*/

    /*std::vector<double> diff_vec;
    mesh.createVector(diff_vec);
    for(unsigned int k=nLocalBegin;k<nLocalEnd;k++)
        diff_vec[k]=dx_funcVal[k]-dx_funcVal1[k];

    mesh.vectorToVTK(diff_vec,"diff");*/

    /*if(rank==0)
        for(unsigned int k=nLocalBegin;k<nLocalEnd;k++)
            std::cout<<"k: "<<k<<" f(k)= "<<funcVal[k]<<" f'(k)= "<<dx_funcVal1[k]<<" f_ap'(k)= "<<dx_funcVal[k]<< " r: "<<(dx_funcVal[k]/dx_funcVal1[k]) <<std::endl;*/

    /*double l2_uzip=normL2(&(*(funcVal.begin()+nLocalBegin)),&(*(dx_funcVal.begin()+nLocalBegin)),(nLocalEnd-nLocalBegin),comm);
    if(!rank)std::cout<<"l2: uzip: "<<l2_uzip<<std::endl;*/

    if(mesh->isActive())
    {

        double l2_x=normL2(&(*(dx_funcVal1.begin()+nLocalBegin)),&(*(dx_funcVal.begin()+nLocalBegin)),(nLocalEnd-nLocalBegin),mesh->getMPICommunicator());
        double l_inf_x= normLInfty(&(*(dx_funcVal1.begin() + nLocalBegin)), &(*(dx_funcVal.begin() + nLocalBegin)),
                                   (nLocalEnd - nLocalBegin), mesh->getMPICommunicator());

        if(!mesh->getMPIRank()) std::cout<<"l2 norm (x derivative): "<<l2_x<<std::endl;
        if(!mesh->getMPIRank()) std::cout<<"l_inf norm (x derivative): "<<l_inf_x<<std::endl;


        double l2_y=normL2(&(*(dy_funcVal1.begin()+nLocalBegin)),&(*(dy_funcVal.begin()+nLocalBegin)),(nLocalEnd-nLocalBegin),mesh->getMPICommunicator());
        double l_inf_y= normLInfty(&(*(dy_funcVal1.begin() + nLocalBegin)), &(*(dy_funcVal.begin() + nLocalBegin)),(nLocalEnd - nLocalBegin), mesh->getMPICommunicator());

        if(!mesh->getMPIRank()) std::cout<<"l2 norm (y derivative): "<<l2_y<<std::endl;
        if(!mesh->getMPIRank()) std::cout<<"l_inf norm (y derivative): "<<l_inf_y<<std::endl;


        double l2_z=normL2(&(*(dz_funcVal1.begin()+nLocalBegin)),&(*(dz_funcVal.begin()+nLocalBegin)),(nLocalEnd-nLocalBegin),mesh->getMPICommunicator());
        double l_inf_z= normLInfty(&(*(dz_funcVal1.begin() + nLocalBegin)), &(*(dz_funcVal.begin() + nLocalBegin)),(nLocalEnd - nLocalBegin), mesh->getMPICommunicator());

        if(!mesh->getMPIRank()) std::cout<<"l2 norm (z derivative): "<<l2_z<<std::endl;
        if(!mesh->getMPIRank()) std::cout<<"l_inf norm (z derivative): "<<l_inf_z<<std::endl;

    }



    /*if(!rank)
        for(unsigned int k=nLocalBegin;k<nLocalEnd;k++)
            std::cout<<" k: "<<k<<" f(k): "<<funcVal[k]<<" dxf(k): "<<dx_funcVal1[k]<<" Sf(k): "<<dx_funcVal[k]<<std::endl;*/

    delete mesh;

    MPI_Barrier(MPI_COMM_WORLD);
    if(!rank) std::cout<<"stencil applied "<<std::endl;

    MPI_Finalize();

    return 0;

}


