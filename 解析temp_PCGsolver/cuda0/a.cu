#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PP 10

static int c[10][10];
static double d[PP][PP];

void MM();

int main()
{
    int i, j;

    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < 10; j++)
        {
            c[i][j] = i * 10 + j;
            d[i][j] = i * 10 + j;
        }
    }
    MM();

    for (i = 0; i < 2; i++)
    printf("%le\n", sin(0));

    printf("%d\n", PP);
    printf("%le\n", d[2][2]);
}

void MM()
{
    d[2][1] = 99;
    printf("%d\n", c[2][1]);
}