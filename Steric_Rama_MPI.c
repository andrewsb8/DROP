//Current compilation command: mpicc Steric_Rama_MPI.c -o steric dihedralRotation.c vectorCalculus.c readProtein.c stericClash.c -lm
//Current run command: mpirun -n 4 steric GGG_COOH_hydrogens_connect.pdb

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

//my includes
#include "readProtein.h"
#include "dihedralRotation.h"
#include "vector_calculus.h"
#include "stericClash.h"

int main(int argc, char *argv[])
{
  //starting positions in Ramachandran space
  double phi = -179;
  double psi = -179;
  double s_phi; //angles passed to worker processes
  double s_psi;

  //start paralleliztion of rotating bonds
  int myrank, nprocs;
  //number of sub grids along each axis
  double phi_subgrids;
  double psi_subgrids;
  //edge sizes along each axis
  double phi_edge;
  double psi_edge;


  //initiate parallelization
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Status status;

  //Master process reads PDB, calculates the dihedral angles, and distributes portions of Rama space to other nodes
  if(myrank == 0)
  {
    //read and print file name from command line
    printf("filename: %s\n", argv[1]);
    char* filename = argv[1];

    //THIS SECTION: Read a pdb, identify bonds, dihedrals, positions, etc
    struct protein prot;
    readPDB(&prot, filename);

    /*//test section for readPDB function to check things read correctly
    printf("%d\n", prot.atoms[2].atom_number);
    printf("%s\n", prot.atoms[2].atom_type);
    printf("%s\n", prot.atoms[2].residue);
    printf("res num: %d\n", prot.atoms[2].residue_number);
    printf("%f %f %f\n", prot.atoms[2].coordinates[0], prot.atoms[2].coordinates[1], prot.atoms[2].coordinates[2]);
    printf("%s\n", prot.atoms[2].atom_name);
    printf("LN: %d\n", prot.number_of_atoms);
    */

    //***NOTE***: want to consolidate readPDBbonds to be an option
    //for readPDB instead of its own method
    readPDBbonds(&prot, filename);
     //more test outputs
    printf("num res: %d\n", prot.number_of_residues);
    printf("num bonds: %d\n", prot.number_of_bonds);
    printf("%d %d\n", prot.bonds[0].bond_atomNumbers[0], prot.bonds[0].bond_atomNumbers[1]);

    //***NOTE***: again, identifyDihedrals should be a part of readPDB
    identifyDihedrals(&prot);

    printf("num dih: %d\n", prot.number_of_dihedrals);
    for(int i=0; i<prot.number_of_dihedrals; i++)
    {
      for(int j=0; j<4; j++)
      {
        printf("%d ", prot.dihedrals[i].dihedral_atomNumbers[j]);
      }
      printf("\n\n");
    }

    //***NOTE***: This block of code should be moved to the readPDB method
    //especially since it is just storing information "already in" the pdb
    for(int i = 0; i < prot.number_of_dihedrals; i++)
    {
      prot.dihedrals[i].dihedral_angle = calculateDihedral(&prot, i);
      printf("%f\n", prot.dihedrals[i].dihedral_angle);
    }
    printf("\n");

    int check = 0;
    char frame[40];
    sprintf(frame, "%s %d", "Frame ", 0);
    writeXYZ(&prot, "trialanine.xyz", frame, 'm', 0, myrank);
    int j;
    for(int i = 1; i <= 180; i++)
    {
      for(j = 1; j <= 180; j++)
      {
        rotateDihedral(&prot, 3, prot.dihedrals[3].dihedral_angle, 2);
        //printf("results: %d %f %f %f\n",i,calculateDihedral(&prot, 0), tmp, calculateDihedral(&prot, 0)-tmp);
        sprintf(frame, "%s %d", "Frame ", i*j);
        writeXYZ(&prot, "trialanine.xyz", frame, 'm', i*j, myrank);
      }
      //checkClashes(&prot);
      //printf("HERE\n");
      //double tmp = calculateDihedral(&prot, 0);
      rotateDihedral(&prot, 0, prot.dihedrals[0].dihedral_angle, 2);
      //printf("results: %d %f %f %f\n",i,calculateDihedral(&prot, 0), tmp, calculateDihedral(&prot, 0)-tmp);
      sprintf(frame, "%s %d", "Frame ", i*j);
      writeXYZ(&prot, "trialanine.xyz", frame, 'm', i*j, myrank);
      //printXYZ(&prot);
    }

    //Next up: want to adjust all of the angles to the following specifications for 4 nodes
    //-179,-179,-179,-179
    //-179,-179,0,-179
    //-179,0,-179,-179
    //-179,0,0,-179

    /*
    //NEXT SECTION: Start distributing protein configurations for analysis

    //if there are more than two nodes utilized, subgrid generation can be generalized
    if(nprocs > 1)
    {
      int i=0;
      int j,k;

      //Need to split the number of grids into phi,psi space
      //grids in each dimension are equal if grid number is even
      if(nprocs % 2 == 0)
      {
        phi_subgrids = nprocs/2;
        psi_subgrids = 2;

        phi_edge = (double)358.0/phi_subgrids;
        psi_edge = (double)358.0/psi_subgrids;

        for(j = 0; j < phi_subgrids; j++) //phi loop
        {
          for(k = 0; k < psi_subgrids; k++) //psi loop
          {
            s_phi = phi + j*phi_edge;
            s_psi = psi + k*psi_edge;

            i++;
            printf("Grid %d: [%f,%f] to [%f,%f]\n",i,s_phi,s_psi,s_phi+phi_edge,s_psi+psi_edge);

            //rotate protein to desired phi, psi combination for grid

            //here will be an if statement that will allow all but the last subgrid to be sent
            //to worker nodes
            //if(j != phi_subgrids-1 && k != psi_subgrids-1){
            //MPI_Send()
            //}
          }
        }


      }
      //if number of grids is odd
      else
      {
        phi_subgrids = ceil((double)nprocs/2);
        psi_subgrids = 2;

        phi_edge = (double)358.0/phi_subgrids;
        psi_edge = (double)358.0/psi_subgrids;

        for(j = 0; j < phi_subgrids-1; j++) //phi loop
        {
          for(k = 0; k < psi_subgrids; k++) //psi loop
          {
            s_phi = phi + j*phi_edge;
            s_psi = psi + k*psi_edge;

            i++;
            printf("Grid %d: [%f,%f] to [%f,%f]\n",i,s_phi,s_psi,s_phi+phi_edge,s_psi+psi_edge);
          }
        }

        s_phi = phi + (phi_subgrids-1)*phi_edge;
        s_psi = psi;

        printf("Grid %d: [%f,%f] to [%f,%f]\n",nprocs,s_phi,s_psi,s_phi+phi_edge,s_psi+(2*psi_edge));

        //rotate protein to desired phi, psi combination for grid

        //here will be an if statement that will allow all but the last subgrid to be sent
        //to worker nodes
        //if(j != phi_subgrids-1 && k != psi_subgrids-1){
        //MPI_Send()
        //}

      }


    }*/

  }
  else
  {
    printf("workers are currently yeeting\n");
  }

  MPI_Finalize();
  return 0;
}
