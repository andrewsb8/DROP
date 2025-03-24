#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#include "../logging/logging.h"
#include "../readProtein/readProtein.h"
#include "../vectorCalculus/vectorCalculus.h"

double calculateDihedral(struct protein *prot, int dihedralNumber)
{
	double angle;
	double angle_in_rads;
	double signedAngle;
	double dot;
	double norm;
	double sign;

	double *vec1 =
		vectorSubtract(prot->atoms
					   [prot->
						dihedrals[dihedralNumber].dihedral_atomNumbers[0] -
						1].coordinates,
					   prot->atoms[prot->
								   dihedrals
								   [dihedralNumber].dihedral_atomNumbers[1]
								   - 1].coordinates);
	double *vec2 =
		vectorSubtract(prot->atoms
					   [prot->
						dihedrals[dihedralNumber].dihedral_atomNumbers[1] -
						1].coordinates,
					   prot->atoms[prot->
								   dihedrals
								   [dihedralNumber].dihedral_atomNumbers[2]
								   - 1].coordinates);
	double *vec3 =
		vectorSubtract(prot->atoms
					   [prot->
						dihedrals[dihedralNumber].dihedral_atomNumbers[2] -
						1].coordinates,
					   prot->atoms[prot->
								   dihedrals
								   [dihedralNumber].dihedral_atomNumbers[3]
								   - 1].coordinates);

	double *cross1 = crossProduct(vec1, vec2);
	double *cross2 = crossProduct(vec2, vec3);

	double *doublecross = crossProduct(cross1, cross2);

	norm = vectorMagnitude(cross1) * vectorMagnitude(cross2);
	dot = dotProduct(cross1, cross2);
	angle_in_rads = acosl(dot / norm);
	angle = (180 / M_PI) * angle_in_rads;

	double ip = dotProduct(vec1, doublecross);
	if (ip < 0) {
		sign = -1;
	} else {
		sign = 1;
	}
	signedAngle = sign * angle;

	free(vec1);
	free(vec2);
	free(vec3);

	free(cross1);
	free(cross2);
	free(doublecross);

	return signedAngle;
}

//find index associated with the relevant dihedral information
int findDihedral(struct protein *prot, int rnum, char *dtype)
{
	int index = -1;
	for (int i = 0; i < prot->number_of_dihedrals; i++) {
		if (rnum == prot->dihedrals[i].dihedral_resNum &&
			strcmp(dtype, *prot->dihedrals[i].dihedral_angType) == 0) {
			index = i;
		}
	}
	return index;
}

int *findDihedrals(struct protein *prot, int rnum, FILE * log)
{
	int *indices = malloc(sizeof(int) * sizeDihedralList);
	for (int i = 0; i < sizeDihedralList; i++) {
		indices[i] = findDihedral(prot, rnum, DihedralList[i]);
		if (indices[i] == -1 && i > 2) {
			char message[120];
			sprintf(message,
					"Warning: Dihedral type %s in residue %d not found.\n",
					DihedralList[i], rnum);
			drop_warning(log, message);
		} else if (indices[i] == -1 && i < 2) {
			char message[120];
			sprintf(message,
					"ERROR: Dihedral type %s in residue %d not found. Backbone dihedrals are required for this calculation.\n",
					DihedralList[i], rnum);
			drop_fatal(log, message);
		}
	}
	return indices;
}

void
updatePositions(struct protein *prot, double newPositions[3],
				int atomNumber)
{
	for (int j = 0; j < 3; j++) {
		prot->atoms[atomNumber].coordinates[j] = newPositions[j];
	}
}

void
rotateDihedral(struct protein *prot, int dihedralNumber,
			   double dihedralAngleChange, bool backbone)
{
	//translate all atoms such that the second atom of the dihedral of interest is at the origin
	int atom_to_origin =
		prot->dihedrals[dihedralNumber].dihedral_atomNumbers[1];
	int atom_rotation_index =
		prot->dihedrals[dihedralNumber].dihedral_atomNumbers[2];

	double translation[3];		//copy of translation vector
	for (int i = 0; i < 3; i++) {
		translation[i] = prot->atoms[atom_to_origin - 1].coordinates[i];
	}

	//translate all atoms
	for (int i = 0; i < prot->number_of_atoms; i++) {
		double *tmp =
			vectorSubtract(translation, prot->atoms[i].coordinates);
		//equate tmp values to prot->atoms[i].coordinates
		updatePositions(prot, tmp, i);
		free(tmp);
	}

	//rotate about z axis to x-z plane
	double angleToXZPlane =
		atanl((prot->atoms[atom_rotation_index - 1].coordinates[1]) /
			  (prot->atoms[atom_rotation_index - 1].coordinates[0]));
	for (int i = 0; i < prot->number_of_atoms; i++) {
		double *tmp = vectorRotate(prot->atoms[i].coordinates, 2,
								   -angleToXZPlane);
		updatePositions(prot, tmp, i);
		free(tmp);
	}

	//rotate about y axis so that the axis between the middle two atoms of the dihedral is aligned with the z axis
	double angleToZAxis =
		atanl((prot->atoms[atom_rotation_index - 1].coordinates[0]) /
			  (prot->atoms[atom_rotation_index - 1].coordinates[2]));
	for (int i = 0; i < prot->number_of_atoms; i++) {
		double *tmp =
			vectorRotate(prot->atoms[i].coordinates, 1, -angleToZAxis);
		updatePositions(prot, tmp, i);
		free(tmp);
	}

	//if dihedral is aligned with the z-axis from 0 to -z rather than 0 to z, reverse rotation direction
	//preserves direction of rotation without needing extra rotation to align with positive z axis
	if (prot->atoms
		[prot->dihedrals[dihedralNumber].dihedral_atomNumbers[2] -
		 1].coordinates[2] < 0) {
		dihedralAngleChange = -dihedralAngleChange;
	}
	//rotate all atoms after bond in question about the z axis by the desired change here
	if (backbone) {
		for (int i = atom_rotation_index - 1; i < prot->number_of_atoms;
			 i++) {
			double *tmp = vectorRotate(prot->atoms[i].coordinates, 2,
									   (M_PI / 180.0) *
									   dihedralAngleChange);
			updatePositions(prot, tmp, i);
			free(tmp);
		}
	} else						//rotate side chain instead
	{
		//bool to find sc atoms in dihedral. don't want to rotate whole side chain, only atoms which come after relevant dihedral
		int found = 0;

		for (int k = 0;
			 k <
			 prot->residues[prot->dihedrals[dihedralNumber].
							dihedral_resNum - 1].num_sc_atoms; k++) {
			//residue/dihedral is identified by third atom in dihedral, consistent with readProtein.c
			if (found == 0 && prot->dihedrals[dihedralNumber].dihedral_atomNumbers[2] ==
				prot->residues[prot->dihedrals
							   [dihedralNumber].dihedral_resNum -
							   1].sidechain_atoms[k]) {
				found = 1;
			}
			if (found == 1) {
				double *tmp =
					vectorRotate(prot->atoms
								 [prot->residues
								  [prot->dihedrals
								   [dihedralNumber].dihedral_resNum -
								   1].sidechain_atoms[k] -
								  1].coordinates, 2,
								 (M_PI / 180.0) * dihedralAngleChange);
				updatePositions(prot, tmp,
								prot->residues[prot->dihedrals
											   [dihedralNumber].dihedral_resNum
											   - 1].sidechain_atoms[k] -
								1);
				free(tmp);
			}
		}
	}

	//undo the first three steps in reverse order
	for (int i = 0; i < prot->number_of_atoms; i++) {
		double *tmp =
			vectorRotate(prot->atoms[i].coordinates, 1, angleToZAxis);
		updatePositions(prot, tmp, i);
		free(tmp);
	}

	for (int i = 0; i < prot->number_of_atoms; i++) {
		double *tmp =
			vectorRotate(prot->atoms[i].coordinates, 2, angleToXZPlane);
		updatePositions(prot, tmp, i);
		free(tmp);
	}

	for (int i = 0; i < prot->number_of_atoms; i++) {
		double *tmp = vectorAdd(translation, prot->atoms[i].coordinates);
		updatePositions(prot, tmp, i);
		free(tmp);
	}
}

//This method will be used to make successive rotations about the same dihedral
//without translating to protein or molecule back to the original coordinates
//between rotations.
//For now, the rotation methods in this file translate protein back to position in
//box. But that is ok for now.
void
rotateDihedral_noTranslate(struct protein *prot, int dihedralNumber,
						   double dihedralAngleChange)
{
	int atom_to_origin =
		prot->dihedrals[dihedralNumber].dihedral_atomNumbers[1] - 1;
	//rotate all atoms about the z axis by the desired change here
	for (int i = atom_to_origin + 1; i < prot->number_of_atoms; i++) {
		double *tmp = vectorRotate(prot->atoms[i].coordinates, 2,
								   (M_PI / 180.0) * dihedralAngleChange);
		updatePositions(prot, tmp, i);
		free(tmp);
	}
}
