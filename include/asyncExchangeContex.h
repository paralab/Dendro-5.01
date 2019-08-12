//
// Created by milinda on 11/19/18.
//

/**
 * @brief Simple class to manage async data transfer in the ODA class.
 * */


#ifndef DENDRO_5_0_UPDATECTX_H
#define DENDRO_5_0_UPDATECTX_H

#include "mpi.h"
#include <vector>

namespace ot {

    class AsyncExchangeContex {

        private :
            /** pointer to the variable which perform the ghost exchange */
            void* m_uiBuffer;

            /** pointer to the send buffer*/
            void* m_uiSendBuf;

            /** pointer to the send buffer*/
            void* m_uiRecvBuf;

            std::vector<MPI_Request*>  m_uiRequests;

        public:
            /**@brief creates an async ghost exchange contex*/
            AsyncExchangeContex(const void* var)
            {
                m_uiBuffer=(void*)var;
                m_uiSendBuf=NULL;
                m_uiRecvBuf=NULL;
                m_uiRequests.clear();
            }

            /**@brief allocates send buffer for ghost exchange*/
            inline void allocateSendBuffer(size_t bytes)
            {
                m_uiSendBuf=malloc(bytes);
            }

            /**@brief allocates recv buffer for ghost exchange*/
            inline void allocateRecvBuffer(size_t bytes)
            {
                m_uiRecvBuf=malloc(bytes);
            }

            /**@brief allocates send buffer for ghost exchange*/
            inline void deAllocateSendBuffer()
            {
                free(m_uiSendBuf);
                m_uiSendBuf=NULL;
            }

            /**@brief allocates recv buffer for ghost exchange*/
            inline void deAllocateRecvBuffer()
            {
                free(m_uiRecvBuf);
                m_uiRecvBuf=NULL;
            }

            inline void* getSendBuffer() { return m_uiSendBuf;}
            inline void* getRecvBuffer() { return m_uiRecvBuf;}

            inline const void* getBuffer() {return m_uiBuffer;}

            inline std::vector<MPI_Request*>& getRequestList(){ return m_uiRequests;}

            bool operator== (AsyncExchangeContex other) const{
                return( m_uiBuffer == other.m_uiBuffer );
            }

            ~AsyncExchangeContex() {

               /* for(unsigned int i=0;i<m_uiRequests.size();i++)
                {
                    delete m_uiRequests[i];
                    m_uiRequests[i]=NULL;
                }

                m_uiRequests.clear();*/

            }

    };

} //end namespace

#endif //DENDRO_5_0_UPDATECTX_H
