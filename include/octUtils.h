//
// Created by milinda on 9/6/16.
//

#ifndef SFCSORTBENCH_OCTUTILS_H
#define SFCSORTBENCH_OCTUTILS_H

#include "sfcSearch.h"
#include "sfcSort.h"

/**
  @brief A collection of simple functions for manipulating octrees.
Examples: Regular Refinements, Linearizing an octree, I/O,
Nearest Common Ancestor, adding positive boundaries, marking hanging nodes
@author Rahul S. Sampath, rahul.sampath@gmail.com
@author Hari Sundar, hsundar@gmail.com
@author Milinda Fernando ,milinda@cs.utah.edu

 @remarks Most of the functions used for the mesh generation. Most of the implementations are based on the previous implementation of dendro version 4.0

*/


#include "TreeNode.h"
#include <vector>
#include <assert.h>
#include "mpi.h"
#include "parUtils.h"
#include <functional>
#include "mathUtils.h"
#include "sfcSort.h"
#include "dendro.h"

#define OCT2BLK_DECOMP_BLK_FILL_RATIO 0.5 // gurantees how fraction of the block covered by regular octants.
#define OCT2BLK_DECOMP_LEV_GAP 0


/**
 * @author Rahul Sampath
 * @author Hari Sundar
 *
 * @breif add the positive boundary octants for the given octree input.
 * */
void addBoundaryNodesType1(std::vector<ot::TreeNode> &in,
                           std::vector<ot::TreeNode>& bdy,
                           unsigned int dim, unsigned int maxDepth);




int refineOctree(const std::vector<ot::TreeNode> & inp,
                 std::vector<ot::TreeNode> &out);

int refineAndPartitionOctree(const std::vector<ot::TreeNode> & inp,
                             std::vector<ot::TreeNode> &out, MPI_Comm comm);


int createRegularOctree(std::vector<ot::TreeNode>& out, unsigned int lev,unsigned int dim, unsigned int maxDepth, MPI_Comm comm);


/**
 * @author Milinda Fernando
 * @brief ensures that the children of the same parent not partioned arcross different partitions.
 * @param [in] in: input octree.
 * @param [in] comm: MPI communicator.
 * @param [out] in : octree with enforced partition constraint.
 *
 * */
void enforceSiblingsAreNotPartitioned(std::vector<ot::TreeNode> & in,MPI_Comm comm);



/**
 * @author Milinda Fernando
 * @brief Computes blucket splitters for a sorted element array.
 * @param [in] pNodes: const ptr to sorted element list
 * @param [in] lev: level of the parent bucket
 * @param [in] maxDepth: max depth of the octree
 * @param [in] rot_id: rotation id of the bucket
 * @param [in] begin: begining of the bucket.
 * @param [in] end : end of the bucket
 * @param [in] splitters: Splitter counts.
 *
 * */
template <typename  T>
void computeSFCBucketSplitters(const T *pNodes, int lev, unsigned int maxDepth,unsigned char rot_id,DendroIntL &begin, DendroIntL &end, DendroIntL *splitters);




template <typename  T>
void computeSFCBucketSplitters(const T *pNodes, int lev, unsigned int maxDepth,unsigned char rot_id,DendroIntL &begin, DendroIntL &end, DendroIntL *splitters)
{
    if ((lev >= maxDepth) || (begin == end)) {
        // Special Case when the considering level exceeds the max depth.

        for (int ii = 0; ii < NUM_CHILDREN; ii++) {
            int index = (rotations[2 * NUM_CHILDREN * rot_id + ii] - '0');
            int nextIndex = 0;
            if (ii == (NUM_CHILDREN-1))
                nextIndex = ii + 1;
            else
                nextIndex = (rotations[2 * NUM_CHILDREN * rot_id + ii + 1] - '0');

            if (ii == 0) {
                splitters[index] = begin;
                splitters[nextIndex] = end;
                continue;
            }
            splitters[nextIndex] = splitters[index];
        }
        //std::cout<<"End return "<<"maxDepth "<<maxDepth<<" Lev: "<<lev<< " Begin "<<begin <<" End "<<end<<std::endl;
        return;

    }

    register unsigned int cnum;
    register unsigned int cnum_prev=0;
    DendroIntL num_elements=0;
    unsigned int rotation=0;
    DendroIntL count[(NUM_CHILDREN+2)]={};
    //unsigned int pMaxDepth=(lev);
    //pMaxDepth--;
    unsigned int mid_bit = maxDepth - lev - 1;
    count[0]=begin;
    for (DendroIntL i=begin; i<end; ++i) {

        /*cnum = (lev < pNodes[i].getLevel())? 1 +(((((pNodes[i].getZ() & (1u << mid_bit)) >> mid_bit) << 2u) |
                                                  (((pNodes[i].getY() & (1u << mid_bit)) >> mid_bit) << 1u) |
                                                  ((pNodes[i].getX() & (1u << mid_bit)) >>
                                                   mid_bit))):0;*/

        cnum = (lev < pNodes[i].getLevel())? 1 +( (((pNodes[i].getZ() >> mid_bit) & 1u) << 2u) | (((pNodes[i].getY() >> mid_bit) & 1u) << 1u) | ((pNodes[i].getX() >>mid_bit) & 1u)):0;
        count[cnum+1]++;


    }

    DendroIntL loc[NUM_CHILDREN+1];
    T unsorted[NUM_CHILDREN+1];
    unsigned int live = 0;

    //if(count[1]>0) std::cout<<"For rank: "<<rank<<" count [1]:  "<<count[1]<<std::endl;

    for (unsigned int ii = 0; ii < NUM_CHILDREN; ii++) {
        int index = (rotations[2 * NUM_CHILDREN * rot_id + ii] - '0');
        int nextIndex = 0;
        if (ii == (NUM_CHILDREN-1))
            nextIndex = ii + 1;
        else
            nextIndex = (rotations[2 * NUM_CHILDREN * rot_id + ii + 1] - '0');

        if (ii == 0) {
            splitters[index] = begin;
            splitters[nextIndex] = splitters[index]+count[1]+ count[(index+2)]; // number of elements which needs to come before the others due to level constraint.

        }else {
            splitters[nextIndex] = splitters[index] + count[(index + 2)];
        }
        //   if(count[1]>0 & !rank) std::cout<<" Spliter B:"<<index <<" "<<splitters[index]<<" Splitters E "<<nextIndex<<" "<<splitters[nextIndex]<<std::endl;

    }
}








#endif //SFCSORTBENCH_OCTUTILS_H
