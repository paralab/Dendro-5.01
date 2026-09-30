//
// Created by milinda on 7/12/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief Contans utility functions to save/load mesh (octree) and variable list
 *defined on the otree.
 */
//

#ifndef SFCSORTBENCH_CHECKPOINT_H
#define SFCSORTBENCH_CHECKPOINT_H

#include <fstream>
#include <iostream>

#include "TreeNode.h"
#include "json.hpp"
#include "mesh.h"
using json = nlohmann::json;
namespace io {

namespace checkpoint {

/**
 *@brief Write an given octree to a file using binary file format.
 *@param[in] fName: file name
 *@param[in] pNodes: pointer to the begin location of the octree
 *@param[in] num: size of pNodes array or number of elements to write.
 * @note than this will write the elements in the order that they have in pNodes
 *array. binary file format: <number of elements><list of elements ...>
 **/

int writeOctToFile(const char* fName, const ot::TreeNode* pNodes,
                   const unsigned int num);

/**
 *@brief Reads an octree speficied by the file name.
 *param[in] fName: input file name
 *param[out] pNodes: TreeNodes read from the specified file.
 **/

int readOctFromFile(const char* fName, std::vector<ot::TreeNode>& pNodes);

/**
 * @breif: writes the variable vector to a file in binary format.
 * @param[in] fName: input file name.
 * @param[in] pMesh: pointer to the mesh which the variable defined on
 * @param [in] vec: pointer to the begin location of the variable vector.
 * binary file format:
 * <vecsize><localBegin><localEnd><values..>
 *
 * */
template <typename T>
int writeVecToFile(const char* fName, const ot::Mesh* pMesh, const T* vec);

/**
 * @breif: writes the variable vectors to a file in binary format.
 * @param[in] fName: input file name.
 * @param[in] pMesh: pointer to the mesh which the variable defined on
 * @param [in] vec: pointer to variable vector
 * @param [in] numVars: number of variables
 * binary file format:
 * <vecsize><localBegin><localEnd><values..>
 *
 * */
template <typename T>
int writeVecToFile(const char* fName, const ot::Mesh* pMesh, const T** vec,
                   const unsigned int numVars);

/**
 * @breif: reads variable vector from the binary file.
 * @param[in] fName: input file name.
 * @param[in] pMesh: pointer to the mesh which the variable defined on
 * @param [out] vec: pointer to the begin location of the variable
 * vector.(assumes memory allocated)
 * */
template <typename T>
int readVecFromFile(const char* fName, const ot::Mesh* pMesh, T* vec);

/**
 * @breif: reads variable vectors from the binary file.
 * @param[in] fName: input file name.
 * @param[in] pMesh: pointer to the mesh which the variable defined on
 * @param [in] numVars: number of variables
 * @param [out] vec: pointer to the begin location of the variable
 * vector.(assumes memory allocated)
 * */
template <typename T>
int readVecFromFile(const char* fName, const ot::Mesh* pMesh, T** vec,
                    const unsigned int numVars);

}  // end of namespace checkpoint.

}  // end of namespace io

// templated implementations
namespace io {
namespace checkpoint {
template <typename T>
int writeVecToFile(const char* fName, const ot::Mesh* pMesh, const T* vec) {
    dendro::logger::debug("Writing vec to file: ", fName);
    unsigned int numNodes    = pMesh->getNumLocalMeshNodes() +
                               pMesh->getNumPreMeshNodes() +
                               pMesh->getNumPostMeshNodes();
    unsigned int nLocalBegin = pMesh->getNodeLocalBegin();
    unsigned int nLocalEnd   = pMesh->getNodeLocalEnd();

    FILE* outfile            = fopen(fName, "w");
    if (outfile == NULL) {
        std::cout << fName << " file open failed " << std::endl;
        return 1;
    }

    const size_t numLocal = pMesh->getNumLocalMeshNodes();
    bool ok = fwrite(&numNodes, sizeof(unsigned int), 1, outfile) == 1 &&
              fwrite(&nLocalBegin, sizeof(unsigned int), 1, outfile) == 1 &&
              fwrite(&nLocalEnd, sizeof(unsigned int), 1, outfile) == 1;
    if (ok && numNodes > 0)
        ok = fwrite((vec + nLocalBegin), sizeof(T), numLocal, outfile) ==
             numLocal;

    // a full disk can surface only at close, when the buffer is flushed
    ok = (fclose(outfile) == 0) && ok;
    if (!ok) {
        std::cout << fName << " file write failed. " << std::endl;
        return 1;
    }
    dendro::logger::debug("Finished writing vec to file: ", fName);
    return 0;
}

template <typename T>
int writeVecToFile(const char* fName, const ot::Mesh* pMesh, const T** vec,
                   const unsigned int numVars) {
    dendro::logger::debug("Writing vec to file: ", fName);
    unsigned int numNodes    = pMesh->getNumLocalMeshNodes() +
                               pMesh->getNumPreMeshNodes() +
                               pMesh->getNumPostMeshNodes();
    unsigned int nLocalBegin = pMesh->getNodeLocalBegin();
    unsigned int nLocalEnd   = pMesh->getNodeLocalEnd();

    FILE* outfile            = fopen(fName, "w");
    if (outfile == NULL) {
        std::cout << fName << " file open failed " << std::endl;
        return 1;
    }

    const size_t numLocal = pMesh->getNumLocalMeshNodes();
    bool ok = fwrite(&numNodes, sizeof(unsigned int), 1, outfile) == 1 &&
              fwrite(&nLocalBegin, sizeof(unsigned int), 1, outfile) == 1 &&
              fwrite(&nLocalEnd, sizeof(unsigned int), 1, outfile) == 1;
    if (ok && numNodes > 0)
        for (unsigned int i = 0; ok && i < numVars; i++)
            ok = fwrite((vec[i] + nLocalBegin), sizeof(T), numLocal, outfile) ==
                 numLocal;

    // a full disk can surface only at close, when the buffer is flushed
    ok = (fclose(outfile) == 0) && ok;
    if (!ok) {
        std::cout << fName << " file write failed. " << std::endl;
        return 1;
    }
    dendro::logger::debug("Finished writing vec to file: ", fName);
    return 0;
}

template <typename T>
int readVecFromFile(const char* fName, const ot::Mesh* pMesh, T* vec) {
    dendro::logger::debug("Reading vec from file: ", fName);
    unsigned int numNodes, nLocalBegin, nLocalEnd;

    FILE* infile = fopen(fName, "r");
    if (infile == NULL) {
        std::cout << fName << " file open failed " << std::endl;
        return 1;
    }
    if (fread(&numNodes, sizeof(unsigned int), 1, infile) != 1 ||
        fread(&nLocalBegin, sizeof(unsigned int), 1, infile) != 1 ||
        fread(&nLocalEnd, sizeof(unsigned int), 1, infile) != 1) {
        std::cout << fName << " file header is truncated. " << std::endl;
        fclose(infile);
        return 1;
    }

    if (numNodes !=
        (pMesh->getNumPreMeshNodes() + pMesh->getNumLocalMeshNodes() +
         pMesh->getNumPostMeshNodes())) {
        std::cout << fName
                  << " file number of total node mismatched with the mesh. "
                  << std::endl;
        fclose(infile);
        return 1;
    }
    if (nLocalBegin != pMesh->getNodeLocalBegin()) {
        std::cout
            << fName
            << " file local node begin location mismatched with the mesh. "
            << std::endl;
        fclose(infile);
        return 1;
    }
    if (nLocalEnd != pMesh->getNodeLocalEnd()) {
        std::cout << fName
                  << " file local node end location mismatched with the mesh. "
                  << std::endl;
        fclose(infile);
        return 1;
    }

    const size_t numLocal = pMesh->getNumLocalMeshNodes();
    if (numNodes > 0 &&
        fread((vec + nLocalBegin), sizeof(T), numLocal, infile) != numLocal) {
        std::cout << fName << " file is truncated. " << std::endl;
        fclose(infile);
        return 1;
    }

    fclose(infile);
    dendro::logger::debug("Finished reading vec from file: ", fName);
    return 0;
}

template <typename T>
int readVecFromFile(const char* fName, const ot::Mesh* pMesh, T** vec,
                    const unsigned int numVars) {
    dendro::logger::debug("Reading vec from file: ", fName);
    unsigned int numNodes, nLocalBegin, nLocalEnd;

    FILE* infile = fopen(fName, "r");
    if (infile == NULL) {
        std::cout << fName << " file open failed " << std::endl;
        return 1;
    }

    if (fread(&numNodes, sizeof(unsigned int), 1, infile) != 1 ||
        fread(&nLocalBegin, sizeof(unsigned int), 1, infile) != 1 ||
        fread(&nLocalEnd, sizeof(unsigned int), 1, infile) != 1) {
        std::cout << fName << " file header is truncated. " << std::endl;
        fclose(infile);
        return 1;
    }

    if (numNodes !=
        (pMesh->getNumPreMeshNodes() + pMesh->getNumLocalMeshNodes() +
         pMesh->getNumPostMeshNodes())) {
        std::cout << fName
                  << " file number of total node mismatched with the mesh. "
                  << std::endl;
        fclose(infile);
        return 1;
    }
    if (nLocalBegin != pMesh->getNodeLocalBegin()) {
        std::cout
            << fName
            << " file local node begin location mismatched with the mesh. "
            << std::endl;
        fclose(infile);
        return 1;
    }
    if (nLocalEnd != pMesh->getNodeLocalEnd()) {
        std::cout << fName
                  << " file local node end location mismatched with the mesh. "
                  << std::endl;
        fclose(infile);
        return 1;
    }

    const size_t numLocal = pMesh->getNumLocalMeshNodes();
    if (numNodes > 0)
        for (unsigned int i = 0; i < numVars; i++)
            if (fread((vec[i] + nLocalBegin), sizeof(T), numLocal, infile) !=
                numLocal) {
                std::cout << fName << " file is truncated. " << std::endl;
                fclose(infile);
                return 1;
            }

    fclose(infile);
    dendro::logger::debug("Finished reading vec from file: ", fName);
    return 0;
}

}  // end of namespace checkpoint

}  // end of namespace io

#endif  // SFCSORTBENCH_CHECKPOINT_H
