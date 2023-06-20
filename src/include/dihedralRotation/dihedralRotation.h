double calculateDihedral(struct protein *prot, int dihedralNumber);
double determineSign(struct protein *prot, int dihedralNumber); //temp to be removed
int findDihedral(struct protein *prot, int rnum, char *dtype);
void updatePositions(struct protein *prot, double newPositions[3], int atomNumber);
double rotateDihedral(struct protein *prot, int dihedralNumber, double dihedralAngle, double dihedralAngleChange, int bb_or_sc, int chi);
double rotateDihedral_noTranslate(struct protein *prot, int dihedralNumber, double dihedralAngle, double dihedralAngleChange);
