//
// Created by milinda on 5/30/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//


#include "oct2vtk.h"

namespace io
{
    namespace vtk
    {



        void oct2vtk(const ot::TreeNode * pNodes, const unsigned int n, const char * fPrefix, MPI_Comm comm)
        {

            int rank,npes;
            MPI_Comm_rank(comm,&rank);
            MPI_Comm_size(comm,&npes);

            char fname[FNAME_LENGTH];
            char str[2048];
            sprintf(fname,"%s_%d_%d.vtu",fPrefix,rank,npes);

            fp=fopen(fname,"w+");
            if(fp==NULL) {
                std::cout << "rank: " << rank << "[IO Error]: Could not open the vtk file. " << std::endl;
                return ;
            }



            unsigned int dim=m_uiDim;
            DendroIntL num_vertices=n*(1u<<dim);//pMesh->getNumLocalMeshElements()*pMesh->getNumNodesPerElement();
            unsigned int num_cells=n;//pMesh->getNumLocalMeshElements();
            //int num_cells_elements = num_cells * NUM_CHILDREN + num_cells;
            unsigned int nPe=(1u<<dim);//pMesh->getNumNodesPerElement();
            unsigned int eleOrder=1;//pMesh->getElementOrder();
            unsigned int sz;

            write_string("<?xml version=\"1.0\"?>\n");
            write_string("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"");
            //write_string(" compressor=\"vtkZLibDataCompressor\"");
            write_string(" byte_order=\"BigEndian\">\n");
            write_string("  <UnstructuredGrid>\n");


            sprintf(str,"    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%d\">\n",(long long) num_vertices,num_cells);
            write_string(str);


            write_string("      <Points>\n");
            if(!useBinary)
                sprintf(str,"        <DataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            else
                sprintf(str,"        <DataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

            write_string(str);
            write_string("          ");

            for(unsigned int ele=0;ele<n;ele++)
            {
                sz=1u<<(m_uiMaxDepth-pNodes[ele].getLevel());
                assert((sz%eleOrder)==0);
                for(unsigned int k=0;k<(eleOrder+1);k++)
                    for(unsigned int j=0;j<(eleOrder+1);j++)
                        for(unsigned int i=0;i<(eleOrder+1);i++)
                        {

                            if(!useBinary)
                            {
                                //sprintf(str,"          %16.8e %16.8e %16.8e\n",(float)(pElements[ele].getX()+i*(sz/eleOrder)),(float)(pElements[ele].getY()+j*(sz/eleOrder)),(float)(pElements[ele].getZ()+k*(sz/eleOrder)));
                                sprintf(str,"          %d %d %d\n",(pNodes[ele].getX()+i*(sz/eleOrder)),(pNodes[ele].getY()+j*(sz/eleOrder)),(pNodes[ele].getZ()+k*(sz/eleOrder)));
                                write_string(str);
                            }else
                            {
                               /* sprintf(str,"%d\0",(pElements[ele].getX()+i*(sz/eleOrder)));
                                write_string(base64::encode(str).c_str());
                                sprintf(str,"%d\0",(pElements[ele].getY()+j*(sz/eleOrder)));
                                write_string(base64::encode(str).c_str());
                                sprintf(str,"%d\0",(pElements[ele].getZ()+k*(sz/eleOrder)));
                                write_string(base64::encode(str).c_str());*/

                                write_int_base64((pNodes[ele].getX()+i*(sz/eleOrder)));
                                write_int_base64((pNodes[ele].getY()+j*(sz/eleOrder)));
                                write_int_base64((pNodes[ele].getZ()+k*(sz/eleOrder)));

                            }



                        }


            }

           write_string("\n");

            write_string("        </DataArray>\n");
            write_string("      </Points>\n");
            write_string("      <Cells>\n");

            if(!useBinary)
                sprintf(str,"        <DataArray type=\"%s\" Name=\"connectivity\""  " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_ASCII);
            else
                sprintf(str,"        <DataArray type=\"%s\" Name=\"connectivity\""  " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_BINARY);
            write_string(str);

            for (unsigned int i = 0 ; i < num_cells ; i++)
            {

                if(!useBinary)
                {
                    sprintf(str,"          %lld %lld %lld %lld %lld %lld %lld %lld\n",
                            (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0),
                            (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder),
                            (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder),
                            (DendroIntL)(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0),
                            (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0),
                            (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder),
                            (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder),
                            (DendroIntL)(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0));
                    write_string(str);

                }else
                {

                    //write_string("          ");

                   /* sprintf(str,"%d\0",i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0);
                    write_string(base64::encode(str).c_str());

                    sprintf(str,"%d\0",i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder);
                    write_string(base64::encode(str).c_str());

                    sprintf(str,"%d\0",i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder);
                    write_string(base64::encode(str).c_str());

                    sprintf(str,"%d\0",i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0);
                    write_string(base64::encode(str).c_str());

                    sprintf(str,"%d\0",i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0);
                    write_string(base64::encode(str).c_str());

                    sprintf(str,"%d\0",i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder);
                    write_string(base64::encode(str).c_str());

                    sprintf(str,"%d\0",i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder);
                    write_string(base64::encode(str).c_str());

                    sprintf(str,"%d\0",i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0);
                    write_string(base64::encode(str).c_str());*/

                    write_int_base64(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0);
                    write_int_base64(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder);
                    write_int_base64(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder);
                    write_int_base64(i * nPe + 0*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0);

                    write_int_base64(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+0);
                    write_int_base64(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+0*(eleOrder+1)+eleOrder);
                    write_int_base64(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+eleOrder);
                    write_int_base64(i * nPe + eleOrder*(eleOrder+1)*(eleOrder+1)+eleOrder*(eleOrder+1)+0);

                }

            }

            write_string("\n");
            write_string("        </DataArray>\n");


            if(!useBinary)
                sprintf(str,"        <DataArray type=\"%s\" Name=\"offsets\"" " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_ASCII);
            else
                sprintf(str,"        <DataArray type=\"%s\" Name=\"offsets\"" " format=\"%s\">\n",DENDRO_NODE_ID_TYPE,DENDRO_FORMAT_BINARY);

            write_string(str);

            if(!useBinary)
            {
                write_string("         ");
                for (unsigned int il = 1,sk = 1; il <= num_cells; ++il, ++sk) {
                    sprintf(str," %lld",(NUM_CHILDREN*il));
                    write_string(str);
                    if (!(sk % 8) && il != num_cells)
                        write_string("\n         ");
                }
                write_string("\n");
            }else
            {
                //write_string("          ");
                for(unsigned int i=1;i<=num_cells;i++)
                {
                   /* sprintf(str,"%lld\0",(long long )((i)*NUM_CHILDREN));
                    write_string(base64::encode(str).c_str());*/
                    write_int_base64(i*(NUM_CHILDREN));
                }


            }

            write_string("\n");
            write_string("        </DataArray>\n");

            if(!useBinary)
                sprintf (str, "        <DataArray type=\"UInt8\" Name=\"types\"" " format=\"%s\">\n", DENDRO_FORMAT_ASCII);
            else
                sprintf (str, "        <DataArray type=\"UInt8\" Name=\"types\"" " format=\"%s\">\n", DENDRO_FORMAT_BINARY);

            write_string(str);

            if(!useBinary)
            {
                write_string("         ");//fprintf (vtufile, "         ");
                for (unsigned int il = 0, sk = 1; il < num_cells; ++il, ++sk) {
                    sprintf (str, " %d",VTK_HEXAHEDRON);
                    write_string(str);
                    if (!(sk % 20) && il != (num_cells - 1))
                        write_string("\n         ");//fprintf (vtufile, "\n         ");
                }
                write_string("\n");
            }else
            {
                //write_string("          ");
                for(unsigned int i=0;i<num_cells;i++)
                {
                    /*sprintf(str,"%d\0",VTK_HEXAHEDRON);
                    write_string(base64::encode(str).c_str());*/
                    write_int_base64(VTK_HEXAHEDRON);
                }


                write_string("\n");

            }


            write_string("        </DataArray>\n");
            write_string("      </Cells>\n");

            // writing cell data
            write_string("      <CellData>\n");

            if(!useBinary)
                sprintf(str,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            else
                sprintf(str,"        <DataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

            write_string(str);
            if(!useBinary)
            {
                write_string ("         ");
                for (unsigned int il = 0; il < num_cells; ++il) {
                    write_int(rank);//fprintf (vtufile, " %d", mpirank);
                }
                write_string("\n");

            }else
            {
                //write_string ("         ");
                for (unsigned int il = 0; il < num_cells; ++il) {
                    write_int_base64(rank);//fprintf (vtufile, " %d", mpirank);
                }
                write_string("\n");
            }


            write_string("        </DataArray>\n");

            if(!useBinary)
                sprintf(str,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
            else
                sprintf(str,"        <DataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\">\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

            write_string(str);
            if(!useBinary)
            {
                write_string ("         ");
                for (unsigned int il = 0; il < num_cells; ++il) {
                    write_int(pNodes[il].getLevel());//fprintf (vtufile, " %d", mpirank);
                }
                write_string("\n");

            }else
            {
                //write_string ("         ");
                for (unsigned int il = 0; il < num_cells; ++il) {
                    write_int_base64(pNodes[il].getLevel());//fprintf (vtufile, " %d", mpirank);
                }
                write_string("\n");
            }


            write_string("        </DataArray>\n");
            write_string("      </CellData>\n");




            write_string("</Piece>\n");
            write_string("</UnstructuredGrid>\n");
            write_string("</VTKFile>\n");

            fclose(fp);
            fp=NULL;

            if(!rank)
            {
                sprintf(fname,"%s.pvtu",fPrefix);

                fp=fopen(fname,"w+");
                if(fp==NULL) {
                    std::cout << "rank: " << rank << "[IO Error]: Could not open the pvtu file. " << std::endl;
                    return ;
                }


                write_string("<?xml version=\"1.0\"?>\n");
                write_string("<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
                if(useBinary) write_string(" compressor=\"vtkZLibDataCompressor\"");
                write_string(" byte_order=\"BigEndian\">\n");

                write_string("  <PUnstructuredGrid GhostLevel=\"0\">\n");



                write_string("    <PPoints>\n");

                if(!useBinary)
                    sprintf(str,"      <PDataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
                else
                    sprintf(str,"      <PDataArray type=\"%s\" Name=\"Position\""" NumberOfComponents=\"3\" format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

                write_string(str);
                write_string("    </PPoints>\n");

                write_string("    <PCellData>\n");

                if(!useBinary)
                    sprintf(str,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
                else
                    sprintf(str,"        <PDataArray type=\"%s\" Name=\"mpi_rank\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

                write_string(str);

                if(!useBinary)
                    sprintf(str,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_ASCII);
                else
                    sprintf(str,"        <PDataArray type=\"%s\" Name=\"cell_level\"" " format=\"%s\"/>\n",DENDRO_NODE_COORD_TYPE,DENDRO_FORMAT_BINARY);

                write_string(str);
                write_string("    </PCellData>\n");



                for(unsigned int proc=0;proc<npes;proc++)
                {
                    sprintf(str,"<Piece Source=\"%s_%d_%d.vtu\"/>\n",fPrefix,proc,npes);
                    write_string(str);
                }
                write_string("</PUnstructuredGrid>\n");
                write_string("</VTKFile>\n");


                fclose(fp);
                fp=NULL;




            }










        }


    } // end of namespace vtk

} // end of namespace io