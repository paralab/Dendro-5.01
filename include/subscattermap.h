/**
 * @file subscattermap.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief sub scatter map for partial MPI communications. This is an abstract version of the full scatter map to perform partial
 * data exchanges.  (Note that this is implemented for block wise scatter map, which will be used in non uniform timestepping )
 * @version 0.1
 * @date 2019-07-21
 * @copyright Copyright (c) 2019, School of Computing, University of Utah. 
 */

#pragma once
#include "mesh.h"


namespace ot
{
    class SubScatterMap
    {

        private:
            /***@brief mesh data structure */
            ot::Mesh* m_uiMesh;

            /** */
            


    };

}// end of namespace ot





