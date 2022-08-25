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

int checkClashes(struct protein *prot);
double getVDWRadii(struct VDW *radii, char *atom_name, char *atom_name_2);
