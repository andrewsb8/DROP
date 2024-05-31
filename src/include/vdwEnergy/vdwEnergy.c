#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "vdwEnergy.h"
#include "../readProtein/readProtein.h"
#include "../vectorCalculus/vectorCalculus.h"

//units: Angstroms
struct VDW_params sigma = {};
//units: kJ/mol
struct VDW_params epsilon = {};

double getParam(struct VDW_params param, char *atom_name)
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

//Lorentz-Berthelot mixing for epsilon VDW or LJ parameter
double mixedEpsilon(char *atom_one_name, char *atom_two_name)
{
  double eps_atom_one = getParam(&epsilon, atom_one_name);
  double eps_atom_two = getParam(&epsilon, atom_two_name);

  return sqrt(eps_atom_one*eps_atom_two);
}

//Lorentz-Berthelot mixing for sigma VDW or LJ parameter
double mixedSigma(char *atom_one_name, char *atom_two_name)
{
  double sig_atom_one = getParam(&sigma, atom_one_name);
  double sig_atom_two = getParam(&sigma, atom_two_name);

  return 0.5*(sig_atom_one+sig_atom_two);
}

//calculate VDW or LJ Energy for protein
double calculateVDWEnergy(struct protein *prot, FILE *log)
{

}
