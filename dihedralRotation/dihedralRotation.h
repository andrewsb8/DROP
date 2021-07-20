double calculateDihedral(struct protein *prot, int dihedralNumber);
void updatePositions(struct protein *prot, double newPositions[3], int atomNumber);
double rotateDihedral(struct protein *prot, int dihedralNumber, double dihedralAngle, double dihedralAngleChange);
double rotateDihedral3(struct protein *prot, int dihedralNumber, double dihedralAngle, double dihedralAngleChange);
