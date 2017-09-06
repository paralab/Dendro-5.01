//
// Created by milinda on 5/30/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief writes a given octree and variable array to vtk file.
*/
//

#ifndef SFCSORTBENCH_OCT2VTK_H
#define SFCSORTBENCH_OCT2VTK_H

#define FNAME_LENGTH 256
#define VTK_HEXAHEDRON 12

#define DENDRO_NODE_COORD_TYPE "UInt32"
#define DENDRO_NODE_ID_TYPE "UInt64"
#define DENDRO_NODE_VAR_INT "UInt32"
#define DENDRO_NODE_VAR_FLOAT "Float32"
#define DENDRO_NODE_VAR_DOUBLE "Float64"

#define DENDRO_FORMAT_ASCII "ascii"
#define DENDRO_FORMAT_BINARY "binary"


#include "TreeNode.h"
#include <iostream>
#include <fstream>
#include "mpi.h"
#include "base.h"
#include "dendro.h"
#include "assert.h"
#include "sfcSort.h"

namespace io
{
    namespace vtk
    {

        static FILE *fp = NULL;
        static int useBinary = 0;
        static int numInColumn = 0;



        static void end_line(void)
        {
            if (!useBinary)
            {
                char str2[8] = "\n";
                fprintf(fp, str2);
                numInColumn = 0;
            }
        }
        static void new_section(void)
        {
            if (numInColumn != 0)
                end_line();
            numInColumn = 0;
        }


        /**
         * @breif Determines if the machine is little-endian.  If so, then, for binary
           data, it will force the data to be big-endian.
         * @note This assumes that all inputs are 4 bytes long.
         * @author: Hank Childs
         * This was taken from the visit website.
         * */
        static void force_big_endian(unsigned char *bytes)
        {
            static int doneTest = 0;
            static int shouldSwap = 0;
            if (!doneTest)
            {
                int tmp1 = 1;
                unsigned char *tmp2 = (unsigned char *) &tmp1;
                if (*tmp2 != 0)
                    shouldSwap = 1;
                doneTest = 1;
            }

            if (shouldSwap & useBinary)
            {
                unsigned char tmp = bytes[0];
                bytes[0] = bytes[3];
                bytes[3] = tmp;
                tmp = bytes[1];
                bytes[1] = bytes[2];
                bytes[2] = tmp;
            }
        }

        /**
         * @breif Writes an integer to the currently open file.
         * @note This assumes that all inputs are 4 bytes long.
         * @author: Hank Childs
         * This was taken from the visit website.
         * */
        static void write_int(int val)
        {
            if (useBinary)
            {
                force_big_endian((unsigned char *) &val);
                fwrite(&val, sizeof(int), 1, fp);
            }
            else
            {
                char str[128];
                sprintf(str, "%d ", val);
                fprintf(fp, str);
                if (((numInColumn++) % 9) == 8)
                {
                    char str2[8] = "\n";
                    fprintf(fp, str2);
                    numInColumn = 0;
                }
            }
        }

        static void write_int_base64(int val)
        {
            /*char str[256];
            force_big_endian((unsigned char *) &val);
            sprintf(str,"%d",val);*/
            char str[sizeof(int)];
            memcpy(str,&val,sizeof(int));
            std::string encoded_64=base64::encode(str);
            char str1[encoded_64.size()+1];
            memcpy(str1,encoded_64.c_str(),sizeof(char)*encoded_64.size());
            str1[encoded_64.size()]='\0';
            force_big_endian((unsigned char *)str1);
            //std::cout<<"encoded_str: "<<encoded_64<<" str1: "<<str1<<std::endl;
            fwrite(str1,sizeof(char),(encoded_64.size()+1),fp);

        }

        /**
         * @breif Writes an float to the currently open file.
         * @note This assumes that all inputs are 4 bytes long.
         * @author: Hank Childs
         * This was taken from the visit website.
         * */
        static void write_float(float val)
        {
            if (useBinary)
            {
                force_big_endian((unsigned char *) &val);
                fwrite(&val, sizeof(float), 1, fp);
            }
            else
            {
                char str[128];
                sprintf(str, "%20.12e ", val);
                fprintf(fp, str);
                if (((numInColumn++) % 9) == 8)
                {
                    end_line();
                }
            }
        }

        /**
         * @breif Writes an double to the currently open file.
         * @note This assumes that all inputs are 4 bytes long.
         * @author: Hank Childs
         * This was taken from the visit website.
         * */
        static void write_double(double val)
        {
            if (useBinary)
            {
                force_big_endian((unsigned char *) &val);
                fwrite(&val, sizeof(double), 1, fp);
            }
            else
            {
                char str[128];
                sprintf(str, "%20.12e ", val);
                fprintf(fp, str);
                if (((numInColumn++) % 9) == 8)
                {
                    end_line();
                }
            }
        }


        /**
         * @breif Writes an string to the currently open file.
         * @note This assumes that all inputs are 4 bytes long.
         * @author: Hank Childs
         * This was taken from the visit website.
         * */
        static void write_string(const char *str)
        {
          fprintf(fp, str);
        }

        /**
         * @breif Writes vtk header to the  open file.
         * @note This assumes that all inputs are 4 bytes long.
         * @author: Hank Childs
         * This was taken from the visit website.
         * */
        static void write_header(void)
        {
            fprintf(fp, "# vtk DataFile Version 2.0\n");
            fprintf(fp, "Dendro-5.0\n");
            if (useBinary)
                fprintf(fp, "BINARY\n");
            else
                fprintf(fp, "ASCII\n");
        }


        /**
        *@breif Writes the given mesh to a binary vtu (in xml format) file.
        * @param [in] pNodes: input nodes
        * @param [in] n: number of elements
        * @param [in] fPrefix: file prefix name
        * @param [in] comm: mpi communicator.
        * */
        void oct2vtk(const ot::TreeNode * pNodes, const unsigned int n, const char * fPrefix, MPI_Comm comm);


    }
}



#endif //SFCSORTBENCH_OCT2VTK_H
