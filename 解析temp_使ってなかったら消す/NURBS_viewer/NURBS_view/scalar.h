/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#ifndef SCALAR_H
#define SCALAR_H

/*
 * Swap integers
 */
void swapIntegers(int *value1, int *value2);

/*
 * Sort two integers
 */
void sortTwoIntegers(int *value1, int *value2);

/*
 * Sort three integers
 */
void sortThreeIntegers(int *value1, int *value2, int *value3);

/*
 * Sort four integers
 */
void sortFourIntegers(int *value1, int *value2, int *value3, int *value4);

/*
 * Sort two integers with function
 */
void sortTwoIntegersWithFunction(int *value1, int *value2,
                                 int (*compare)(const void *a, const void *b));

/*
 * Sort three integers with function
 */
void sortThreeIntegersWithFunction(int *value1, int *value2, int *value3,
                                   int (*compare)(const void *a, const void *b));

/*
 * Sort four integers with function
 */
void sortFourIntegersWithFunction(int *value1, int *value2, int *value3, int *value4,
                                  int (*compare)(const void *a, const void *b));

#endif  /* SCALAR_H */
