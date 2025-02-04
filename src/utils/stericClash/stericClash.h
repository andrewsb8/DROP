#ifndef STERICCLASH_H_
#define STERICCLASH_H_

struct stericRadii {
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

int countClashes(struct protein *prot, FILE * log, bool list_clashes);
double getRadius(struct stericRadii *radii, char *atom_name,
				 char *atom_name_2);

#endif
