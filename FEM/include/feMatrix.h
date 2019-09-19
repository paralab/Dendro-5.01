//
// Created by milinda on 10/30/18.
//
/**
 * @brief class that derived from abstract class feMat
 * LHS computation of the weak formulation
 * */
#ifndef DENDRO_5_0_FEMATRIX_H
#define DENDRO_5_0_FEMATRIX_H

#include "feMat.h"

template<typename T>
void fem_eMatTogMat(ot::MatRecord *eMat, const T *Imat_p2c, unsigned int pOrder,unsigned int dof=1) {
    const unsigned int n = (pOrder + 1);
    const unsigned int nPe = n * n * n;
    const unsigned int N = nPe * dof;
    const unsigned int dof2=dof*dof;

    // todo : need to fix this properly, with efficient matrix matrix multiplication. 
    T * IT_Ke_I  =new T[N * N];
    
    for(unsigned int di =0; di< dof; di++)
    for(unsigned int dj=0; dj< dof; dj++)
    {
        const unsigned int offset = (di*dof + dj)*nPe*nPe;
        for (unsigned int i = 0; i < nPe; i++)
            for(unsigned int j=0; j < nPe; j++)
            {
                T val = (T)0;
                for(unsigned int k=0;k< nPe; k++)
                {
                    val  += (eMat[offset + i*nPe +k].getMatVal() * Imat_p2c[k * nPe + j]);
                }
                IT_Ke_I[offset + i*nPe +j] = val;
            }
    }
    
    for(unsigned int di =0; di< dof; di++)
    for(unsigned int dj=0; dj< dof; dj++)
    {
        const unsigned int offset = (di*dof + dj)*nPe*nPe;

         for (unsigned int i = 0; i < nPe; i++)
            for(unsigned int j=0; j < nPe; j++)
            {
                eMat[offset + i*nPe +j].setMatValue(IT_Ke_I[offset + i*nPe +j]);
            }
    }
    

    for(unsigned int di =0; di< dof; di++)
    for(unsigned int dj =0; dj< dof; dj++)
    {
         const unsigned int offset = (di*dof + dj)*nPe*nPe;
         for (unsigned int i = 0; i < nPe; i++)
            for(unsigned int j=0; j < nPe; j++)
            {
                T val = (T)0;
                for(unsigned int k=0;k< nPe; k++)
                {
                    val  += ( Imat_p2c[k * nPe + i] * eMat[offset + k*nPe +j].getMatVal());
                }
                IT_Ke_I[offset + i*nPe + j] = val;
            }
    }
    
    for(unsigned int di =0; di< dof; di++)
    for(unsigned int dj=0; dj< dof; dj++)
    {
         const unsigned int offset = (di*dof + dj)*nPe*nPe;
         for (unsigned int i = 0; i < nPe; i++)
            for(unsigned int j=0; j < nPe; j++)
            {
                eMat[offset + i*nPe + j].setMatValue(IT_Ke_I[offset + i*nPe + j]);
            }
    }


    delete [] IT_Ke_I; 
}


template<typename T>
class feMatrix : public feMat {

protected:

    /**@brief number of dof*/
    unsigned int m_uiDof;

    /**@brief element nodal vec in */
    VECType *m_uiEleVecIn;

    /***@brief element nodal vecOut */
    VECType *m_uiEleVecOut;

    /** elemental coordinates */
    double *m_uiEleCoords;

public:
    /**
     * @brief constructs an FEM stiffness matrix class.
     * @param[in] da: octree DA
     * */
    feMatrix(ot::DA *da, unsigned int dof = 1);

    ~feMatrix();

    /**@brief Computes the LHS of the weak formulation, normally the stifness matrix times a given vector.
      * @param [in] in input vector u
      * @param [out] out output vector Ku
      * @param [in] default parameter scale vector by scale*Ku
    * */
    virtual void matVec(const VECType *in, VECType *out, double scale = 1.0);


    /**@brief Computes the elemental matvec
      * @param [in] in input vector u
      * @param [out] out output vector Ku
      * @param [in] default parameter scale vector by scale*Ku
    **/
    virtual void elementalMatVec(const VECType *in, VECType *out, double *coords = NULL, double scale = 1.0) = 0;


#ifdef BUILD_WITH_PETSC

    /**@brief Computes the LHS of the weak formulation, normally the stifness matrix times a given vector.
      * @param [in] in input vector u
      * @param [out] out output vector Ku
      * @param [in] default parameter scale vector by scale*Ku
    * */
    virtual void matVec(const Vec &in, Vec &out, double scale = 1.0);

    /**
     * @brief Performs the matrix assembly.
     * @param [in/out] J: Matrix assembled
     * @param [in] mtype: Matrix type
     * when the function returns, J is set to assembled matrix
     **/
    virtual bool getAssembledMatrix(Mat *J, MatType mtype);


#endif

    /**@brief static cast to the leaf node of the inheritance*/
    T &asLeaf() { return static_cast<T &>(*this); }


    /**
     * @brief executed just before  the matVec loop in matvec function
     * @param[in] in : input Vector
     * @param[out] out: output vector
     * @param[in] scale: scalaing factror
     **/

    bool preMatVec(const VECType *in, VECType *out, double scale = 1.0) {
        return asLeaf().preMatVec(in, out, scale);
    }


    /**@brief executed just after the matVec loop in matvec function
     * @param[in] in : input Vector
     * @param[out] out: output vector
     * @param[in] scale: scalaing factror
     * */

    bool postMatVec(const VECType *in, VECType *out, double scale = 1.0) {
        return asLeaf().postMatVec(in, out, scale);
    }

    /**@brief executed before the matrix assembly */
    bool preMat() {
        return asLeaf().preMat();
    }

    /**@brief executed after the matrix assembly */
    bool postMat() {
        return asLeaf().postMat();
    }

    /**
     * @brief Compute the elemental Matrix.
     * @param[in] eleID: element ID
     * @param[in] coords : elemental coordinates
     * @param[out] records: records corresponding to the elemental matrix.
     *
     * */
    void getElementalMatrix(unsigned int eleID, std::vector<ot::MatRecord> &records, double *coords) {
        return asLeaf().getElementalMatrix(eleID, records, coords);
    }


};

template<typename T>
feMatrix<T>::feMatrix(ot::DA *da, unsigned int dof) : feMat(da) {
    m_uiDof = dof;
    const unsigned int nPe = m_uiOctDA->getNumNodesPerElement();
    m_uiEleVecIn = new VECType[m_uiDof * nPe];
    m_uiEleVecOut = new VECType[m_uiDof * nPe];

    m_uiEleCoords = new double[m_uiDim * nPe];

}

template<typename T>
feMatrix<T>::~feMatrix() {
    delete[] m_uiEleVecIn;
    delete[] m_uiEleVecOut;
    delete[] m_uiEleCoords;

    m_uiEleVecIn = NULL;
    m_uiEleVecOut = NULL;
    m_uiEleCoords = NULL;

}

template<typename T>
void feMatrix<T>::matVec(const VECType *in, VECType *out, double scale) {

    VECType *_in = NULL;
    VECType *_out = NULL;

    if (!(m_uiOctDA->isActive()))
        return;


    preMatVec(in, out, scale);

    m_uiOctDA->nodalVecToGhostedNodal(in, _in, false, m_uiDof);
    m_uiOctDA->createVector(_out, false, true, m_uiDof);

    VECType *val = new VECType[m_uiDof];
    for (unsigned int var = 0; var < m_uiDof; var++)
        val[var] = (VECType) 0;

    m_uiOctDA->setVectorByScalar(_out, val, false, true, m_uiDof);

    delete[] val;

    const ot::Mesh * pMesh = m_uiOctDA->getMesh();
    const unsigned int nPe = pMesh->getNumNodesPerElement();
    double * qMat = new double[nPe*nPe];
    const unsigned int * e2n_cg = &(*(pMesh->getE2NMapping().begin()));
    const ot::TreeNode* allElements = &(*(pMesh->getAllElements().begin()));
    
    m_uiOctDA->readFromGhostBegin(_in, m_uiDof);
    

    const unsigned int totalNodalSize = m_uiOctDA->getTotalNodalSz();
    for (m_uiOctDA->init<ot::DA_FLAGS::INDEPENDENT>(); m_uiOctDA->curr() < m_uiOctDA->end<ot::DA_FLAGS::INDEPENDENT>(); m_uiOctDA->next<ot::DA_FLAGS::INDEPENDENT>()) {

        m_uiOctDA->getElementNodalValues(_in, m_uiEleVecIn, m_uiOctDA->curr(), m_uiDof);
        const unsigned int currentId = m_uiOctDA->curr();
        pMesh->getElementQMat(m_uiOctDA->curr(),qMat,true);
        m_uiOctDA->getElementalCoords(m_uiOctDA->curr(), m_uiEleCoords);
        elementalMatVec(m_uiEleVecIn, m_uiEleVecOut, m_uiEleCoords, scale);

        for (unsigned int dof = 0; dof < m_uiDof; dof++) {
          for (unsigned int i = 0; i < nPe; i++) {
            VECType ss = (VECType) 0;
            for (unsigned int j = 0; j < nPe; j++) {
              ss += qMat[j * nPe + i] * m_uiEleVecOut[dof*nPe + j];
            }
            _out[dof*totalNodalSize + e2n_cg[currentId * nPe + i]] += (VECType) ss;

          }
        }

    }

    m_uiOctDA->readFromGhostEnd(_in, m_uiDof);
    
    const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd = pMesh -> getElementLocalEnd();

    for (m_uiOctDA->init<ot::DA_FLAGS::W_DEPENDENT>(); m_uiOctDA->curr() < m_uiOctDA->end<ot::DA_FLAGS::W_DEPENDENT>(); m_uiOctDA->next<ot::DA_FLAGS::W_DEPENDENT>()) {

        // temporary fix to skip ghost writable.     
        if( m_uiOctDA->curr()< eleLocalBegin || m_uiOctDA->curr()>=eleLocalEnd )
            continue;

        m_uiOctDA->getElementNodalValues(_in, m_uiEleVecIn, m_uiOctDA->curr(), m_uiDof);
        const unsigned int currentId = m_uiOctDA->curr();
        pMesh->getElementQMat(m_uiOctDA->curr(),qMat,true);
        m_uiOctDA->getElementalCoords(m_uiOctDA->curr(), m_uiEleCoords);
        elementalMatVec(m_uiEleVecIn, m_uiEleVecOut, m_uiEleCoords, scale);

        for (unsigned int dof = 0; dof < m_uiDof; dof++) {
          for (unsigned int i = 0; i < nPe; i++) {
            VECType ss = (VECType) 0;
            for (unsigned int j = 0; j < nPe; j++) {
              ss += qMat[j * nPe + i] * m_uiEleVecOut[dof*nPe + j];
            }
            _out[dof*totalNodalSize + e2n_cg[currentId * nPe + i]] += (VECType) ss;

          }
        }

    }


    delete [] qMat;

    // accumilate from write ghost. 
    m_uiOctDA->writeToGhostsBegin(_out,m_uiDof);
    m_uiOctDA->writeToGhostsEnd(_out,ot::DA_FLAGS::WriteMode::ADD_VALUES,m_uiDof);
    
    m_uiOctDA->ghostedNodalToNodalVec(_out, out, true, m_uiDof);

    m_uiOctDA->destroyVector(_in);
    m_uiOctDA->destroyVector(_out);

    postMatVec(in, out, scale);


    return;

}

#ifdef BUILD_WITH_PETSC

template<typename T>
void feMatrix<T>::matVec(const Vec &in, Vec &out, double scale) {

    const PetscScalar *inArry = NULL;
    PetscScalar *outArry = NULL;

    VecGetArrayRead(in, &inArry);
    VecGetArray(out, &outArry);

    matVec(inArry, outArry, scale);

    VecRestoreArrayRead(in, &inArry);
    VecRestoreArray(out, &outArry);

}


template<typename T>
bool feMatrix<T>::getAssembledMatrix(Mat *J, MatType mtype) {
    

    if(m_uiOctDA->isActive())
    {
        MatZeroEntries(*J);
        std::vector<ot::MatRecord> records;
        //
        const unsigned int eleOrder = m_uiOctDA->getElementOrder();
        const unsigned int npe_1d = eleOrder + 1;
        const unsigned int npe_2d = (eleOrder + 1) * (eleOrder + 1);
        const unsigned int nPe = (eleOrder + 1) * (eleOrder + 1) * (eleOrder + 1);
        
        const ot::Mesh *pMesh = m_uiOctDA->getMesh();
        const ot::TreeNode *allElements = &(*(pMesh->getAllElements().begin()));
        
        DendroScalar *p2cEleMat = new DendroScalar[nPe * nPe];
        
        preMat();
        unsigned int nCount = 0;

        double *coords = new double[m_uiDim * nPe];
        
        bool faceHang[NUM_FACES];
        bool edgeHang[NUM_EDGES];
        unsigned int cnumFace[NUM_FACES];
        unsigned int cnumEdge[NUM_EDGES];

        const unsigned int * e2n_cg = &(*(pMesh->getE2NMapping().begin()));
        const unsigned int * e2n_dg = &(*(pMesh->getE2NMapping_DG().begin()));
        const unsigned int * e2e = &(*(pMesh->getE2EMapping().begin()));

        const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
        const unsigned int eleLocalEnd = pMesh -> getElementLocalEnd();

        for (m_uiOctDA->init<ot::DA_FLAGS::LOCAL_ELEMENTS>();m_uiOctDA->curr() < m_uiOctDA->end<ot::DA_FLAGS::LOCAL_ELEMENTS>(); m_uiOctDA->next<ot::DA_FLAGS::LOCAL_ELEMENTS>()) {
            
            const unsigned int currentId = m_uiOctDA->curr();
            m_uiOctDA->getElementalCoords(currentId, coords);
            getElementalMatrix(m_uiOctDA->curr(), records, coords);
            const unsigned int cnum = allElements[currentId].getMortonIndex();
            
            pMesh->getElementQMat(m_uiOctDA->curr(),p2cEleMat,true);
            //printArray_2D(p2cEleMat,nPe,nPe);


            fem_eMatTogMat(&(*(records.begin() + nCount)), p2cEleMat, eleOrder,m_uiDof);
            nCount += (m_uiDof*nPe*nPe*m_uiDof);
            if (records.size() > 500) {
                m_uiOctDA->petscSetValuesInMatrix(*J, records, m_uiDof, ADD_VALUES);
                nCount = 0;
            }

        }//end writable

        m_uiOctDA->petscSetValuesInMatrix(*J, records, m_uiDof, ADD_VALUES);

        postMat();

        MatAssemblyBegin(*J, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(*J, MAT_FINAL_ASSEMBLY);

        delete[] p2cEleMat;
        delete[] coords;

        PetscFunctionReturn(0);
    }

    
}


#endif


#endif //DENDRO_5_0_FEMATRIX_H
