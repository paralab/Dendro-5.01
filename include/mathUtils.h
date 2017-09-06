//
// Created by milinda on 3/31/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains useful math helper routines
*/
//

#ifndef SFCSORTBENCH_MATHUTILS_H
#define SFCSORTBENCH_MATHUTILS_H

#include <iostream>
#include <cmath>
#include <math.h>
#include "parUtils.h"

/**
 * @brief computes the l2 norm between two vectors.
 * @param[in] vec1 input vector 1
 * @param[in] vec2 input vector 2
 * @param[in] dimension of the vector.
 * @param[in] comm communicator
 * @return l2 norm between two vectors.
 * */
template <typename T>
T normL2(T * vec1,T* vec2, unsigned int n,MPI_Comm comm);

/**
 * @brief computes the l2 norm between two vectors.
 * @param[in] vec1 input vector 1
 * @param[in] vec2 input vector 2
 * @param[in] dimension of the vector.
 * @param[in] comm communicator
 * @return l1 norm between two vectors.
 * */
template <typename T>
T normLInfty(T *vec1, T *vec2, unsigned int n, MPI_Comm comm);


/**
 * @brief computes the l2 norm between two vectors.
 * @param[in] vec1 input vector 1
 * @param[in] vec2 input vector 2
 * @param[in] dimension of the vector.
 * @return l1 norm between two vectors.
 * */
template <typename T>
T normLInfty(T *vec1, T *vec2, unsigned int n);



/**
 * @brief computes the l2 norm of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @param[in] comm communicator
 * @return l2 norm of vec.
 * */
template <typename T>
T normL2(T * vec,unsigned int n, MPI_Comm comm);


/**
 * @brief computes the l2 norm between two vectors.
 * @param[in] vec1 input vector 1
 * @param[in] vec2 input vector 2
 * @param[in] dimension of the vector.
 * @return l2 norm between two vectors.
 * */
template <typename T>
T normL2(T * vec1,T* vec2, unsigned int n);

/**
 * @brief computes the l2 norm of a vector.
 * @param[in] vec input vector
 * @param[in] dimension of the vector.
 * @return l2 norm of vec.
 * */
template <typename T>
T normL2(T * vec,unsigned int n);






#include "mathUtils.tcc"


#endif //SFCSORTBENCH_MATHUTILS_H
