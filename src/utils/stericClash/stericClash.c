#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

//keep imports in this order to avoid type conflict error
#include "../readProtein/readProtein.h"
#include "stericClash.h"
#include "../vectorCalculus/vectorCalculus.h"

//may turn different radii into an option to be selected

//normally allowed via Ramachandran JMB 1963 in Angstroms
//Values:           C-C, C-O, C-N, C-H, O-O, O-N, O-H, N-N, N-H, H-H
//struct VDW radii = {3.2, 2.8, 2.9, 2.4, 2.8, 2.7, 2.4, 2.7, 2.4, 2.0};
//inner limit allowed via Ramachandran JMB 1963 in Angstroms
struct stericRadii radii =
	{ 3.0, 2.7, 2.8, 2.2, 2.7, 2.6, 2.2, 2.6, 2.2, 1.9 };

int countClashes(struct protein *prot, FILE * log, bool list_clashes)
{
	double distance;
	int count = 0;

	if (list_clashes) {
		fprintf(log,
				"Clash: (AtomType AtomNumber)-(AtomType AtomNumber) [Allowed Distance, Distance in Structure]\n");
	}
	//loop through all atom combinations i and j where j > i
	for (int i = 0; i < prot->number_of_atoms; i++) {
		for (int j = i + 1; j < prot->number_of_atoms; j++) {
			//check that at least 4 number of covalent bonds are between the atoms being compared
			if (prot->atoms[i].covalent_bondArray[j - i - 1] > 3) {
				//if this condition is satisfied, check distance between atoms, and the difference of van der waals distances
				double *bond_vector =
					vectorSubtract(prot->atoms[i].coordinates,
								   prot->atoms[j].coordinates);
				distance = vectorMagnitude(bond_vector);
				free(bond_vector);
				double min_distance_allowed = getRadius(&radii,
														prot->atoms
														[i].atom_name,
														prot->atoms
														[j].atom_name);
				if (distance < min_distance_allowed) {
					if (list_clashes) {
						fprintf(log,
								"Clash: (%s%d)-(%s%d) [%f, %f]\n",
								prot->atoms[i].atom_name,
								prot->atoms[i].atom_number,
								prot->atoms[j].atom_name,
								prot->atoms[j].atom_number,
								min_distance_allowed, distance);
					}
					count += 1;
				}
			}
		}
	}

	return count;
}

double getRadius(struct stericRadii *radii, char *atom_name,
				 char *atom_name_2)
{
	switch (*atom_name) {
	case 'N':
		switch (*atom_name_2) {
		case 'N':
			return radii->NN;
		case 'H':
			return radii->NH;
		case 'C':
			return radii->CN;
		case 'O':
			return radii->ON;
		}
	case 'C':
		switch (*atom_name_2) {
		case 'N':
			return radii->CN;
		case 'H':
			return radii->CH;
		case 'C':
			return radii->CC;
		case 'O':
			return radii->CO;
		}
	case 'O':
		switch (*atom_name_2) {
		case 'N':
			return radii->ON;
		case 'H':
			return radii->OH;
		case 'C':
			return radii->CO;
		case 'O':
			return radii->OO;
		}
	case 'H':
		switch (*atom_name_2) {
		case 'N':
			return radii->NH;
		case 'H':
			return radii->HH;
		case 'C':
			return radii->CH;
		case 'O':
			return radii->OH;
		}
	default:
		printf("No atom name defined for %s.\n", atom_name);
		break;
	}
}
