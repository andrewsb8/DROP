#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

//keep imports in this order to avoid type conflict error
#include "../readProtein/readProtein.h"
#include "vdwEnergy.h"
#include "../vectorCalculus/vectorCalculus.h"
#include "../stericClash/stericClash.h"

//atom order: C, N, O, H
//atom types from CHARMM36m: carbonyl C (peptide backbone), amide Nitrogen, carbonyl oxygen, nonpolar H
//units: Angstroms
struct VDW_params sigma = {3.56359, 3.29632, 3.02905, 2.35197};
//units: kJ/mol
struct VDW_params epsilon = {0.46024, 0.83680, 0.50208, 0.09204};

//copied from stericClash.c, inner radii for steric clashes
//struct VDW radii2 = {3.0, 2.7, 2.8, 2.2, 2.7, 2.6, 2.2, 2.6, 2.2, 1.9};

double getParam(struct VDW_params *param, char *atom_name)
{
  switch(*atom_name)
  {
    case 'N' :
      return param->N;
    case 'C' :
      return param->C;
    case 'O' :
      return param->O;
    case 'H' :
      return param->H;
    default :
      printf("No atom name defined for %s.\n", *atom_name);
      break;
  }
}

//Lorentz-Berthelot mixing for epsilon VDW/LJ parameter
double mixedEpsilon(char *atom_one_name, char *atom_two_name)
{
  double eps_atom_one = getParam(&epsilon, atom_one_name);
  double eps_atom_two = getParam(&epsilon, atom_two_name);

  return sqrt(eps_atom_one*eps_atom_two);
}

//Lorentz-Berthelot mixing for sigma VDW/LJ parameter
double mixedSigma(char *atom_one_name, char *atom_two_name)
{
  double sig_atom_one = getParam(&sigma, atom_one_name);
  double sig_atom_two = getParam(&sigma, atom_two_name);

  return 0.5*(sig_atom_one+sig_atom_two);
}

double pairVDWEnergy(double distance, double epsilon, double sigma)
{
  return 4*epsilon*( pow(sigma/distance, 12) - pow(sigma/distance, 6) );
}

//calculate VDW or LJ Energy for protein
double calculateVDWEnergy(struct protein *prot, double gamma, FILE *log)
{
  double distance;
  double mixed_sigma;
  double mixed_epsilon;
  double energy = 0;

  //loop through all atom combinations i and j where j > i
  for(int i = 0; i < prot->number_of_atoms; i++)
  {
    for(int j = i+1; j < prot->number_of_atoms; j++)
    {
      //check that at least 4 number of covalent bonds are between the atoms being compared
      if(prot->atoms[i].covalent_bondArray[j-i-1] > 3)
      {
        double *bond_vector = vectorSubtract(prot->atoms[i].coordinates,prot->atoms[j].coordinates);
        distance = vectorMagnitude(bond_vector);
        free(bond_vector);
        if(distance < gamma * getVDWRadii(&radii, prot->atoms[i].atom_name, prot->atoms[j].atom_name) && ((isBackbone(prot->atoms[i].atom_type) && !isBackbone(prot->atoms[j].atom_type)) || (!isBackbone(prot->atoms[i].atom_type) && isBackbone(prot->atoms[j].atom_type)) ))
        {
          return NAN;
        }
        mixed_sigma = mixedSigma(prot->atoms[i].atom_name, prot->atoms[j].atom_name);
        mixed_epsilon = mixedEpsilon(prot->atoms[i].atom_name, prot->atoms[j].atom_name);
        energy += pairVDWEnergy(distance, mixed_epsilon, mixed_sigma);
      }
    }
  }
  return energy;
}
