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

int countClashes(struct protein *prot);
double getVDWRadii(struct VDW *radii, char *atom_name, char *atom_name_2);

#endif
