//
// Created by milinda on 9/6/16.
//
/**
  @brief A collection of simple functions for manipulating octrees.
Examples: Regular Refinements, Linearizing an octree, I/O,
Nearest Common Ancestor, adding positive boundaries, marking hanging nodes
@author Rahul S. Sampath, rahul.sampath@gmail.com
@author Hari Sundar, hsundar@gmail.com
@author Milinda Fernando ,milinda@cs.utah.edu

 @remarks Most of the functions used for the mesh generation. Most of the implementations are based on the previous implementation of dendro version 4.0

*/

#include "octUtils.h"


// This will add boundary nodes and will also embed the octree one level higher
// to enable the addition of the boundary nodes. The positive boundary nodes
// are also marked as BOUNDARY.
void addBoundaryNodesType1(std::vector<ot::TreeNode> &in,
                           std::vector<ot::TreeNode>& bdy,
                           unsigned int dim, unsigned int maxDepth) {

    assert(bdy.empty());

    for (unsigned int i = 0; i < in.size(); i++) {
        // get basic info ...
        unsigned int d   = in[i].getLevel();
        unsigned int x = in[i].getX();
        unsigned int y = in[i].getY();
        unsigned int z = in[i].getZ();

        unsigned char bdyFlags;
        // check if this is a boundary octant or not ...
        if ( in[i].isBoundaryOctant(ot::TreeNode::POSITIVE, &bdyFlags) ) {
            // bdy flags tells us which octants to add ...

            //NOTE: == is important since a&(b+c) will be true if
            //a=b, a=c and a=b+c

            // +x and more ... add additional octant in +x dir
            if ( bdyFlags & ot::TreeNode::X_POS_BDY ) {
                bdy.push_back(ot::TreeNode( (1u << maxDepth), y, z, (d+1) |
                                                                    ot::TreeNode::BOUNDARY, dim, maxDepth+1));
            }

            // +y and more ... add additional octant in +y dir
            if ( bdyFlags & ot::TreeNode::Y_POS_BDY ) {
                bdy.push_back(ot::TreeNode(x, (1u << maxDepth), z,
                                           (d+1) | ot::TreeNode::BOUNDARY, dim, maxDepth+1));
            }

            // +z and more ... add additional octant in +z dir
            if ( bdyFlags & ot::TreeNode::Z_POS_BDY ) {
                bdy.push_back(ot::TreeNode(x, y, (1u << maxDepth),
                                           (d+1) | ot::TreeNode::BOUNDARY, dim, maxDepth+1));
            }

            //+x+y and more
            if ( (bdyFlags & (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY)) ==
                 (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY) ) {
                bdy.push_back(ot::TreeNode((1u << maxDepth),(1u << maxDepth), z,
                                           (d+1) | ot::TreeNode::BOUNDARY, dim, maxDepth+1));
            }

            //+x+z and more
            if ( (bdyFlags & (ot::TreeNode::X_POS_BDY + ot::TreeNode::Z_POS_BDY)) ==
                 (ot::TreeNode::X_POS_BDY + ot::TreeNode::Z_POS_BDY) ) {
                bdy.push_back(ot::TreeNode((1u << maxDepth), y, (1u << maxDepth),
                                           (d+1) | ot::TreeNode::BOUNDARY, dim, maxDepth+1));
            }

            //+y+z and more
            if ( (bdyFlags & (ot::TreeNode::Y_POS_BDY + ot::TreeNode::Z_POS_BDY)) ==
                 (ot::TreeNode::Y_POS_BDY + ot::TreeNode::Z_POS_BDY) ) {
                bdy.push_back(ot::TreeNode(x, (1u << maxDepth),(1u << maxDepth),
                                           (d+1) | ot::TreeNode::BOUNDARY, dim, maxDepth+1));
            }

            // if global corner ...
            //+x +y and +z only
            if ( (bdyFlags & (ot::TreeNode::X_POS_BDY + ot::TreeNode::Y_POS_BDY +
                              ot::TreeNode::Z_POS_BDY)) == (ot::TreeNode::X_POS_BDY +
                                                            ot::TreeNode::Y_POS_BDY +  ot::TreeNode::Z_POS_BDY) ) {
                bdy.push_back(ot::TreeNode((1u << maxDepth), (1u << maxDepth), (1u << maxDepth),
                                           (d+1) | ot::TreeNode::BOUNDARY, dim, maxDepth+1));
            }
        }//end if boundary

        // Embed the actual octant in one level higher ...
        in[i] = ot::TreeNode(x, y, z, d+1, dim, maxDepth+1);

    }//end for i

    // A Parallel Sort for the bdy nodes follows in the constructor.
    //Then in and bdy will be merged.

}//end function


int refineOctree(const std::vector<ot::TreeNode> & inp,
                 std::vector<ot::TreeNode> &out) {
    out.clear();
    for(unsigned int i = 0; i < inp.size(); i++) {
        if(inp[i].getLevel() < inp[i].getMaxDepth()) {
            inp[i].addChildren(out);
        } else {
            out.push_back(inp[i]);
        }
    }
    return 1;
}//end function


int refineAndPartitionOctree(const std::vector<ot::TreeNode> & inp,
                             std::vector<ot::TreeNode> &out, MPI_Comm comm) {
    refineOctree(inp,out);
    par::partitionW<ot::TreeNode>(out, NULL,comm);
    return 1;
}//end function

int createRegularOctree(std::vector<ot::TreeNode>& out, unsigned int lev,
                        unsigned int dim, unsigned int maxDepth, MPI_Comm comm) {
    ot::TreeNode root(0,0,0,0,dim,maxDepth);
    out.clear();
    int rank;
    MPI_Comm_rank(comm,&rank);
    if(!rank) {
        out.push_back(root);
    }
    for(int i = 0; i < lev; i++) {
        std::vector<ot::TreeNode> tmp;
        refineAndPartitionOctree(out,tmp,comm);
        out = tmp;
        tmp.clear();
    }
    return 1;
}



void enforceSiblingsAreNotPartitioned(std::vector<ot::TreeNode> & in,MPI_Comm comm)
{

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    unsigned int prev,next;
    (rank<(npes-1)) ? next=rank+1: next=0;
    (rank>0) ? prev=rank-1: prev=(npes-1);

    if(in.size()<NUM_CHILDREN) return ;

    unsigned int blCount=0;
    int sendCount=0;
    int recvCount=0;

    MPI_Request req;
    MPI_Status status;

    ot::TreeNode maxNode=in.back();
    ot::TreeNode prev_max;
    if(rank<(npes-1)) sendCount=1;
    if(rank>0) recvCount=1;

    par::Mpi_Sendrecv(&maxNode,sendCount,next,0,&prev_max,recvCount,prev,0,comm,&status);

    sendCount=0;
    recvCount=0;

    for(unsigned int elem=0;elem<NUM_CHILDREN;elem++) {
        if ((elem < in.size()) && (in.front().getParent() == in[elem].getParent()))
            sendCount++;
        else
            break;
    }

    if((rank>0) && (sendCount==1) && (in.front().getParent()!=prev_max.getParent())) sendCount=0;
    if(sendCount==NUM_CHILDREN || rank==0) sendCount=0;

    //std::cout<<"m_uiRank: "<<m_uiRank<<" sendCount: "<<sendCount<<std::endl;
    par::Mpi_Sendrecv(&sendCount,1,prev,1,&recvCount,1,next,1,comm,&status);
    /*
    MPI_Isend(&sendCount,1,MPI_INT,prev,0,m_uiComm,&req);
    MPI_Recv(&recvCount,1,MPI_INT,next,0,m_uiComm,&status);*/

    //std::cout<<"m_uiRank: "<<m_uiRank<<" recvCount: "<<recvCount<<std::endl;

    assert(in.size()>=sendCount);
    std::vector<ot::TreeNode> recvBuffer;
    recvBuffer.resize(recvCount);
    par::Mpi_Sendrecv(&(*(in.begin())),sendCount,prev,2,&(*(recvBuffer.begin())),recvCount,next,2,comm,&status);
    //std::cout<<"rank: "<<m_uiRank<<" send recv ended "<<std::endl;

    if(sendCount)
    {
        std::vector<ot::TreeNode> tmpElements;
        std::swap(in,tmpElements);
        in.clear();
        in.resize(tmpElements.size()-sendCount);
        for(unsigned int ele=sendCount;ele<tmpElements.size();ele++)
            in[ele-sendCount]=tmpElements[ele];

        tmpElements.clear();
    }

    for(unsigned int ele=0;ele<recvBuffer.size();ele++)
        in.push_back(recvBuffer[ele]);

    recvBuffer.clear();

    assert(seq::test::isUniqueAndSorted(in));

    return ;

}