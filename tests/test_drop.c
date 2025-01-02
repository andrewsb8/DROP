#include <check.h>
#include <stdlib.h>
#include <stdio.h>

START_TEST(test_version_output) {
    FILE *fp;
    char path[1035];

    fp = popen("./drop --version 2>&1", "r");
    if (fp == NULL) {
        ck_abort_msg("Failed to run drop");
    }

    fgets(path, sizeof(path)-1, fp);
    pclose(fp);
    ck_assert_str_eq(path, "DROP Version 2025.1-dev\n");
}
END_TEST

Suite *drop_test_suite(void)
{
    Suite *s;
    TCase *tc_core;

    s = suite_create("drop_binary_tests");
    tc_core = tcase_create("Core");

    tcase_add_test(tc_core, test_version_output);
    suite_add_tcase(s, tc_core);

    return s;
}

int main(void)
{
    int num_failed;
    Suite *s;
    SRunner *sr;

    s = drop_test_suite();
    sr = srunner_create(s);

    srunner_run_all(sr, CK_NORMAL);
    num_failed = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (num_failed == 0) ? 0 : 1;
}
