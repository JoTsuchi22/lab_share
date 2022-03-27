#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <time.h>

static int c[100][2];

void MM()
{
    printf("%d\n", c[2][1]);
}

int main()
{
    // int i = 0;
    // printf("%d\n", i);
    for (int k = 0; k < 100; k++)
    {
        c[k][1] = k;
        // printf("%d\n", c[k][1]);
    }
    MM();
}