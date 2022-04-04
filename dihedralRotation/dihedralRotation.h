double calculateDihedral(struct protein *prot, int dihedralNumber);
void updatePositions(struct protein *prot, double newPositions[3], int atomNumber);
double rotateDihedral(struct protein *prot, int dihedralNumber, double dihedralAngle, double dihedralAngleChange, int bb_or_sc);
double rotateDihedral_noTranslate(struct protein *prot, int dihedralNumber, double dihedralAngle, double dihedralAngleChange);
