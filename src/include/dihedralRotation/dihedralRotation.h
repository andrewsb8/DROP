#ifndef DIHEDRALROTATION_H_
#define DIHEDRALROTATION_H_

double calculateDihedral(struct protein *prot, int dihedralNumber);
double determineSign(struct protein *prot, int dihedralNumber); //temp to be removed
int findDihedral(struct protein *prot, int rnum, char *dtype, FILE *log);
void updatePositions(struct protein *prot, double newPositions[3], int atomNumber);
double rotateDihedral(struct protein *prot, int dihedralNumber, double dihedralAngleChange, bool backbone);
double rotateDihedral_noTranslate(struct protein *prot, int dihedralNumber, double dihedralAngleChange);

#endif
