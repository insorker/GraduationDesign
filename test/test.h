#ifndef TEST_H
#define TEST_H

#include <stdio.h>

#define PRINT_TEST_FILE() printf("==== TEST FILE: %s ====\n", __FILE__);
#define PRINT_TEST_FUNC() printf("----- TEST FUNCTION: %s -----\n", __func__);

#endif
