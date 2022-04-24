#include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <math.h>
// #include <assert.h>
// #include <time.h>


#include "myheader.h"
extern double ZETAWIN[K][K];

void test(int i)
{
	int j;
	printf("%d\n", i);
	printf("%d\n", K);
	for (j = 0; j < K; j++)
	{
		for (int k = 0; k < K; k++)
		{
			ZETAWIN[j][k] = j * K + k;
			if (j == 3)
			printf("ZETAWIN[%d][%d] = %le\n", j, k, ZETAWIN[j][k]);
		}
	}
}