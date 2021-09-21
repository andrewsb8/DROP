struct VDW
{
  double H;
  double C;
  double N;
  double O;
};

int checkClashes(struct protein *prot);
double getVDWRadii(struct VDW *radii, char *atom_name);
