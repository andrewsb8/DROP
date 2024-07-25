#ifndef STERICCLASH_H_
#define STERICCLASH_H_

struct VDW
{
  double CC;
  double CO;
  double CN;
  double CH;

  double OO;
  double ON;
  double OH;

  double NN;
  double NH;

  double HH;
};

extern struct VDW radii;

int countClashes (struct protein *prot, FILE * log, bool list_clashes);
double getVDWRadii (struct VDW *radii, char *atom_name, char *atom_name_2);

#endif
