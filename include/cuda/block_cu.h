//
// Created by milinda on 8/22/18.
//
/**
 * @author Milinda Fernando
 * School of Computing, University of Utah
 * @brief simplified ot::Block class for cuda.
 * */

#ifndef DENDRO_5_0_BLOCK_CU_H
#define DENDRO_5_0_BLOCK_CU_H

#include "dendro.h"

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif


namespace cuda
{

        class _Block
        {

            private:

                    /**min point of the block*/
                    double m_uiPtMin[3];

                    /**max point of the block*/
                    double m_uiPtMax[3];

                    /**offset for unzip vars*/
                    unsigned int m_uiOffset;

                    /**blundary flag*/
                    unsigned int m_uiBFlag;

                    /**block size in 3D*/
                    unsigned int m_uiSz[3];

                    /**align the block dimensions*/
                    unsigned int m_uiAlignSz[3];

                    /**dx size*/
                    double m_uiDx[3];


            public:

            /**default constructor*/
            CUDA_CALLABLE_MEMBER _Block()
            {

                m_uiPtMin[0]=0.0;
                m_uiPtMin[1]=0.0;
                m_uiPtMin[2]=0.0;


                m_uiPtMax[0]=0.0;
                m_uiPtMax[1]=0.0;
                m_uiPtMax[2]=0.0;

                m_uiOffset=0;

                m_uiBFlag=0;

                m_uiSz[0]=0;
                m_uiSz[1]=0;
                m_uiSz[2]=0;

                m_uiDx[0]=0.0;
                m_uiDx[1]=0.0;
                m_uiDx[2]=0.0;

                m_uiAlignSz[0]=0;
                m_uiAlignSz[1]=0;
                m_uiAlignSz[2]=0;




            }

            /**creates a valid dendro block on the gpu. */
            CUDA_CALLABLE_MEMBER _Block(const double * p_ptmin, const double * p_ptmax, unsigned int p_offset, unsigned int p_bflag, const unsigned int * p_sz, const double * p_dx )
            {

                m_uiPtMin[0]=p_ptmin[0];
                m_uiPtMin[1]=p_ptmin[1];
                m_uiPtMin[2]=p_ptmin[2];


                m_uiPtMax[0]=p_ptmax[0];
                m_uiPtMax[1]=p_ptmax[1];
                m_uiPtMax[2]=p_ptmin[2];

                m_uiOffset=p_offset;

                m_uiBFlag=p_bflag;

                m_uiSz[0]=p_sz[0];
                m_uiSz[1]=p_sz[1];
                m_uiSz[2]=p_sz[2];

                m_uiDx[0]=p_dx[0];
                m_uiDx[1]=p_dx[1];
                m_uiDx[2]=p_dx[2];

                ((m_uiSz[0] & ((1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG)-1))==0)? m_uiAlignSz[0]=m_uiSz[0] : m_uiAlignSz[0]=((m_uiSz[0]/(1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG))+1)*(1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG);
                m_uiAlignSz[1]=m_uiSz[1];
                m_uiAlignSz[2]=m_uiSz[2];
                //((m_uiSz[1] & ((1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG)-1))==0)? m_uiAlignSz[1]=m_uiSz[1] : m_uiAlignSz[1]=((m_uiSz[1]/(1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG))+1)*(1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG);
                //((m_uiSz[2] & ((1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG)-1))==0)? m_uiAlignSz[2]=m_uiSz[2] : m_uiAlignSz[2]=((m_uiSz[2]/(1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG))+1)*(1u<<DENDRO_BLOCK_ALIGN_FACTOR_LOG);



            }

            /**@returns pt min of the block used in bssn computations*/
            CUDA_CALLABLE_MEMBER const double * getPtMin() const
            {
                return m_uiPtMin;
            }


            /**@returns pt max of the block used in bssn computations*/
            CUDA_CALLABLE_MEMBER const double * getPtMax() const
            {
                return m_uiPtMax;
            }

            /**@returns get offset*/
            CUDA_CALLABLE_MEMBER unsigned int getOffset() const
            {
                return m_uiOffset;
            }


            /**@returns get bflag*/
            CUDA_CALLABLE_MEMBER unsigned int getBFlag() const
            {
                return m_uiBFlag;
            }

            /**@returns get spatial dx*/
            CUDA_CALLABLE_MEMBER const unsigned int * getSz() const
            {
                return  m_uiSz;
            }


            /**@returns get spatial dx*/
            CUDA_CALLABLE_MEMBER const unsigned int * getAlignedSz() const
            {
                return  m_uiAlignSz;
            }


            /**@returns get aligned size of the block*/
            CUDA_CALLABLE_MEMBER const unsigned int getAlignedBlockSz() const
            {
                return  m_uiAlignSz[0]*m_uiAlignSz[1]*m_uiAlignSz[2];
            }


            /**@returns get blocks sz*/

            CUDA_CALLABLE_MEMBER const double * getDx() const
            {
                return m_uiDx;
            }


        };


} // end of namespace cuda























#endif //DENDRO_5_0_BLOCK_CU_H



