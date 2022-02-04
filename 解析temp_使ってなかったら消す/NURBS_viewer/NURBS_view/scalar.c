/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#include "scalar.h"

static int compareIntegers(const void *a, const void *b);

/*
 * Swap integers
 */
void swapIntegers(int *value1, int *value2)
{
    int temp = *value1;
    *value1  = *value2;
    *value2  = temp;
}

/*
 * Sort two integers
 */
void sortTwoIntegers(int *value1, int *value2)
{
    sortTwoIntegersWithFunction(value1, value2,
                                compareIntegers);
}

/*
 * Sort three integers
 */
void sortThreeIntegers(int *value1, int *value2, int *value3)
{
    sortThreeIntegersWithFunction(value1, value2, value3,
                                  compareIntegers);
}

/*
 * Sort four integers
 */
void sortFourIntegers(int *value1, int *value2, int *value3, int *value4)
{
    sortFourIntegersWithFunction(value1, value2, value3, value4,
                                 compareIntegers);
}

/*
 * Sort two integers with function
 */
void sortTwoIntegersWithFunction(int *value1, int *value2,
                                 int (*compare)(const void *a, const void *b))
{
    if (compare(value1, value2) > 0)
        swapIntegers(value1, value2);
}

/*
 * Sort three integers with function
 */
void sortThreeIntegersWithFunction(int *value1, int *value2, int *value3,
                                   int (*compare)(const void *a, const void *b))
{
    if (compare(value1, value2) > 0)
        swapIntegers(value1, value2);
    if (compare(value2, value3) > 0)
        swapIntegers(value2, value3);

    if (compare(value1, value2) > 0)
        swapIntegers(value1, value2);
}

/*
 * Sort four integers with function
 */
void sortFourIntegersWithFunction(int *value1, int *value2, int *value3, int *value4,
                                  int (*compare)(const void *a, const void *b))
{
    if (compare(value1, value2) > 0)
        swapIntegers(value1, value2);
    if (compare(value2, value3) > 0)
        swapIntegers(value2, value3);
    if (compare(value3, value4) > 0)
        swapIntegers(value3, value4);

    if (compare(value1, value2) > 0)
        swapIntegers(value1, value2);
    if (compare(value2, value3) > 0)
        swapIntegers(value2, value3);

    if (compare(value1, value2) > 0)
        swapIntegers(value1, value2);
}

/*
 * Compare integers
 */
static int compareIntegers(const void *a, const void *b)
{
    return *(const int *)a - *(const int *)b;
}
