#include <stdio.h>

#include "myheader.h"
extern double ZETAWIN[K][K];

void test()
{
	int j;
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