#include <check.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

#include "../../src/utils/readProtein/readProtein.h"

// test to validate parsing of pdb ATOM line including:
// - highest atom and residue number
// - four letter atom types
// - four character residues (special or nonstandard residues)
// - two letter atom names at end
// - the chain identifier is ignored
START_TEST(test_pdb_atom_parse) {
    char *line = "ATOM  99999 xCNH GLYNX9999      19.593  19.795  21.155  1.00  0.00           CN";
    struct protein prot;
    size_t size = sizeof(struct _atoms);
	prot.atoms = (struct _atoms *) malloc(size);
    readPDBAtom(&prot, line, 0);
    ck_assert_int_eq(prot.atoms[0].atom_number, 99999);
    ck_assert_str_eq(prot.atoms[0].atom_type, "xCNH");
    ck_assert_str_eq(prot.atoms[0].residue, "GLYN");
    ck_assert_int_eq(prot.atoms[0].residue_number, 9999);
    ck_assert_float_eq(prot.atoms[0].coordinates[0], 19.593);
    ck_assert_float_eq(prot.atoms[0].coordinates[1], 19.795);
    ck_assert_float_eq(prot.atoms[0].coordinates[2], 21.155);
    ck_assert_str_eq(prot.atoms[0].atom_name, "CN");
}
END_TEST

Suite *pdb_io(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("pdb_io");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_pdb_atom_parse);
    suite_add_tcase(s, tc_core);

    return s;
}
