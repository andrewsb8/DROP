#ifndef VDWENERGY_H_
#define VDWENERGY_H_

struct VDW_params
{
  double C;
  double N;
  double O;
  double H;
};

double getParam(struct VDW_params param, char * atom_name);
double mixedEpsilon(char * atom_one_name, char * atom_two_name);
double mixedSigma(char * atom_one_name, char * atom_two_name);
double calculateVDWEnergy(struct protein *prot, FILE *log);

#endif
