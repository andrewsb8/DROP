//Current compilation command: mpicc -g Steric_SC_Rama_MPI.c -o steric dihedralRotation/dihedralRotation.c vectorCalculus/vectorCalculus.c readProtein/readProtein.c stericClash/stericClash.c -lm
//Current run command: mpirun -n 4 steric /path/GGG_COOH_hydrogens_connect.pdb

/* This executable is going to scan the top left of Ramachandran space
(pPII and beta) of a central amino acid (Alanine) of a tripeptide. It will
consider steric interactions and count the number of accessible side chain states.
The idea is to calculate a Ramachandran landscape that is essentially a
Potential of Mean Force (PMF) landscape based on accessible states of the side
chain in context of the backbone configuration.

Author: Brian Andrews
Institution: Drexel University
Last Date Modified: 4/4/2022
*/

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

//my includes
#include "readProtein/readProtein.h"
#include "dihedralRotation/dihedralRotation.h"
#include "vectorCalculus/vectorCalculus.h"
#include "stericClash/stericClash.h"

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

    printf("num atoms: %d\n", prot.number_of_atoms);
    printf("num res: %d\n", prot.number_of_residues);
    printf("num bonds: %d\n", prot.number_of_bonds);
    printf("%d %d\n", prot.bonds[0].bond_atomNumbers[0], prot.bonds[0].bond_atomNumbers[1]);

    printf("num dih: %d\n", prot.number_of_dihedrals);
    for(int i=0; i<prot.number_of_dihedrals; i++)
    {
      for(int j=0; j<4; j++)
      {
        printf("%d ", prot.dihedrals[i].dihedral_atomNumbers[j]);
      }
      printf("\n\n");
    }

    //get phi and psi of central residue to phi = -179 and psi = 179
    //may move this section to an individual method
    //need to have a better idea of how to identify dihedrals from arrays of atom numbers
    double tmp = calculateDihedral(&prot, 3);
    printf("%f\n", tmp);
    rotateDihedral(&prot, 3, prot.dihedrals[3].dihedral_angle, 2, 1, 0);
    double tmp2 = calculateDihedral(&prot, 3);
    printf("%f\n", tmp2);
    if(tmp2 - tmp > 0)
    {
      rotateDihedral(&prot, 3, prot.dihedrals[3].dihedral_angle, (179-tmp2), 1, 0);
    }
    else
    {
      rotateDihedral(&prot, 3, prot.dihedrals[3].dihedral_angle, tmp2-179, 1, 0);
    }

    printf("%f\n", calculateDihedral(&prot, 3));

    //phi
    tmp = calculateDihedral(&prot, 0);
    printf("%f\n", tmp);
    rotateDihedral(&prot, 0, prot.dihedrals[0].dihedral_angle, 2, 1, 0);
    tmp2 = calculateDihedral(&prot, 0);
    printf("%f\n", tmp2);
    if(tmp2 - tmp > 0)
    {
      rotateDihedral(&prot, 0, prot.dihedrals[0].dihedral_angle, 179-tmp2-2, 1, 0);
    }
    else
    {
      rotateDihedral(&prot, 0, prot.dihedrals[0].dihedral_angle, tmp2-179, 1, 0);
    }

    printf("%f\n", calculateDihedral(&prot, 0));

    int check = 0;
    int clash;
    int allowed = 0;
    char frame[40];
    int j;
    sprintf(frame, "%s %d", "Frame ", 0);
    writeXYZ(&prot, "trialanine_ILE2_BBandSC.xyz", frame, 'm', 0, myrank);
    //FILE *free_spaces;
    for(int k = 1; k <= 1; k++) //phi -179 to -41 in 2 degree intervals (69)
    {
      for(int j = 1; j <= 2; j++) //psi 179 to 99 (41)
      {
        for(int h = 1; h <= 180; h++) //chi 2 all of space (180)
        {
          for(int i = 1; i <= 180; i++) //chi 2 all of space
          {
            rotateDihedral(&prot, 5, prot.dihedrals[5].dihedral_angle, 2, 0, 2);
            sprintf(frame, "%s %d", "Frame ", i);
            writeXYZ(&prot, "trialanine_ILE2_BBandSC.xyz", frame, 'm', i, myrank);
            //if(checkClashes(&prot) == 0)
            //{
            //  allowed += 1;  //increment number of allowed states
            //}
            //printf("%f %d\n", calculateDihedral(&prot, 5), checkClashes(&prot));
            //printXYZ(&prot);
          }
          rotateDihedral(&prot, 5, prot.dihedrals[5].dihedral_angle, 2, 0, 2);
          rotateDihedral(&prot, 4, prot.dihedrals[5].dihedral_angle, 2, 0, 1);
        }
        printf("%f %f %d\n", calculateDihedral(&prot, 0), calculateDihedral(&prot, 3), allowed);
        rotateDihedral(&prot, 4, prot.dihedrals[5].dihedral_angle, 2, 0, 1);
        rotateDihedral(&prot, 3, prot.dihedrals[3].dihedral_angle, -2, 1, 0);

        allowed = 0; //reset allowed states
      }

      rotateDihedral(&prot, 3, prot.dihedrals[3].dihedral_angle, 82, 1, 0); //reset psi to 179
      rotateDihedral(&prot, 0, prot.dihedrals[0].dihedral_angle, 2, 1, 0);

    }

    /*
    //NEXT SECTION: Start distributing protein configurations for analysis

    //if there are more than two nodes utilized, subgrid generation can be generalized
    if(nprocs > 1)
    {


    }*/

  }
  else
  {
    printf("workers are currently yeeting\n");
  }

  MPI_Finalize();
  return 0;
}
