/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains useful math helper routines
*
*/



template <typename T>
T normL2(T * vec1,T* vec2, unsigned int n,MPI_Comm comm)
{

    int rank,npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    /*unsigned int n_global;
    par::Mpi_Reduce(&n,&n_global,1,MPI_SUM,0,comm);*/

    T l2=0.0;
    for(unsigned int i=0;i<n;i++)
        l2+=pow((vec1[i]-vec2[i]),2);

    T l2_sum=0;
    par::Mpi_Reduce(&l2,&l2_sum,1,MPI_SUM,0,comm);

    return sqrt(l2_sum);

}

template <typename T>
T normLInfty(T *vec1, T *vec2, unsigned int n, MPI_Comm comm)
{
    int rank,npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    assert(n>0);
    /*unsigned int n_global;
    par::Mpi_Reduce(&n,&n_global,1,MPI_SUM,0,comm);*/

    T l1=fabs(vec1[0]-vec2[0]);
    for(unsigned int i=1;i<n;i++)
        if(l1<(fabs(vec1[i]-vec2[i]))) l1=fabs(vec1[i]-vec2[i]);

    T l1_max=0;
    par::Mpi_Reduce(&l1,&l1_max,1,MPI_MAX,0,comm);

    return (l1_max);
}

template <typename T>
T normLInfty(T *vec1, T *vec2, unsigned int n)
{
 assert(n>0);
    /*unsigned int n_global;
    par::Mpi_Reduce(&n,&n_global,1,MPI_SUM,0,comm);*/

    T l1=fabs(vec1[0]-vec2[0]);
    for(unsigned int i=1;i<n;i++)
        if(l1<(fabs(vec1[i]-vec2[i]))) l1=fabs(vec1[i]-vec2[i]);

    return l1;

}


template <typename T>
T normL2(T * vec,unsigned int n,MPI_Comm comm)
{
    T l2=0;
    for(unsigned int i=0;i<n;i++)
        l2+=pow((vec[i]),2);

    T l2_sum=0;
    par::Mpi_Reduce(&l2,&l2_sum,1,MPI_SUM,0,comm);

    return sqrt(l2_sum);

}


template <typename T>
T normL2(T * vec1,T* vec2, unsigned int n)
{



    T l2=0.0;
    for(unsigned int i=0;i<n;i++)
        l2+=pow((vec1[i]-vec2[i]),2);

    return sqrt(l2);

}

template <typename T>
T normL2(T * vec,unsigned int n)
{
    T l2=0;
    for(unsigned int i=0;i<n;i++)
        l2+=pow((vec[i]),2);


    return sqrt(l2);

}