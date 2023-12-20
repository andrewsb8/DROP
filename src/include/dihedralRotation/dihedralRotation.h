#include <stdbool.h>

double calculateDihedral(struct protein *prot, int dihedralNumber);
int findDihedral(struct protein *prot, int rnum, char *dtype);
void updatePositions(struct protein *prot, double newPositions[3], int atomNumber);
double rotateDihedral(struct protein *prot, int dihedralNumber, double dihedralAngleChange, bool backbone);
double rotateDihedral_noTranslate(struct protein *prot, int dihedralNumber, double dihedralAngleChange);
