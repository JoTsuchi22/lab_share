#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static double AA[10][10];
static double L[10][10];
static double d[10];
static int max_A = 0;

void Cholesky(int n)
{
    L[0][0] = AA[0][0];
    d[0] = 1.0 / L[0][0];
 
    for(int i = 1; i < n; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            if(fabs(AA[i][j]) < 1.0e-10) continue;
 
            double lld = AA[i][j];
            for(int k = 0; k < j; ++k)
            {
                lld -= L[i][k]*L[j][k]*d[k];
            }
            L[i][j] = lld;
        }
 
        d[i] = 1.0 / L[i][i];
    }
}


void Cholesky_v2(int n, double *Aarray, int *rowarray, int *colarray, double *LTarray, double *ddarray)
{
    // L と d ではなく L^T と d を求める関数であることに注意(上三角)
    // ↑ 元のK_whole_val (double *A) に入れるため   (不完全コレスキー分解)
    LTarray[0] = Aarray[0];
    dd[0] = 1.0 / LTarray[0];
    int i, j, k;
    int ddcount = 1;
    int icount = 1;
    for (i = 1; i < max_A; i++)
    {
        for (j = 0; j < col[icount]; j++)
        {
            double lld = Aarray[icount]; //= A[j][i]のこと j:row[icount], i:col[icount]


            for (k = 0; k < row[icount]; k++)
            {
                if ()
                {
                    lld -= LTarray[] * LTarray[] * ddarray[];
                }
            }
            
            LTarray[icount] = lld;

            if (row[icount] == col[icount])
            {
                ddarray[ddcount] = 1.0 / LTarray[icount];
                ddcount++;
            }

            icount++;
        }
    }
 
    for(int i = 1; i < n; ++i)
    {
        for(int j = 0; j <= i; ++j)
        {
            // if(fabs(AA[i][j]) < 1.0e-10) continue;
 
            double lld = AA[i][j];
            for(int k = 0; k < j; ++k)
            {
                lld -= L[i][k]*L[j][k]*d[k];
            }
            L[i][j] = lld;
        }
 
        d[i] = 1.0 / L[i][i];
    }
}


int main()
{
    int ndof = 10;
    double *A = (double *)malloc(sizeof(double) * ndof * ndof);
    int *row = (int *)malloc(sizeof(int) * ndof * ndof);
    int *col = (int *)malloc(sizeof(int) * ndof * ndof);
    double *LT = (double *)malloc(sizeof(double) * ndof * ndof);
    double *dd = (double *)malloc(sizeof(double) * ndof);

    int i, j, k;

    int count = 1;
    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < 10; j++)
        {
            if (i <= j)
            {
                AA[i][i] = 6.0;
                if (i + 1 < 10)
                AA[i][i + 1] = -1.0;
                if (i + 3 < 10)
                AA[i][i + 3] = -1.0;
                count++;
            }
            else
            {
                AA[i][j] = AA[j][i];
            }
        }
    }

    int icount = 0;
    for (i = 0; i < 10; i++)
    {
        int temp = 0;
        A[icount] = 6.0;
        col[icount] = i;
        row[icount] = i;
        if (i + 1 < 10)
        {
            A[icount+1] = -1.0;
            col[icount+1] = i + 1;
            row[icount+1] = i;
            temp++;
        }
        if (i + 3 < 10)
        {
            A[icount+2] = -1.0;
            col[icount+2] = i + 3;
            row[icount+2] = i;
            temp++;
        }
        icount += temp + 1;
    }

    max_A = icount;

    icount = 0;
    printf("A\n");
    printf("A[icount]\trow[icount]\tcol[icount]\n");
    for (i = 0; i < max_A; i++)
    {
        printf("%le\t%d\t%d\n", A[icount], row[icount], col[icount]);
        icount++;
    }
    printf("\n");

    printf("AA\n");
    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < 10; j++)
        {
            printf("%le ", AA[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    Cholesky(ndof);

    printf("L\n");
    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < 10; j++)
        {
            printf("%le ", L[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("d\n");
    for (i = 0; i < 10; i++)
    {
        printf("%le ", d[i]);
    }
    printf("\n\n");

    double BB1[10][10], BB2[10][10];
    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < 10; j++)
        {
            BB1[i][j] = 0.0;
            BB2[i][j] = 0.0;
        }
    }

    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < 10; j++)
        {
            for (k = 0; k < 10; k++)
            {
                if (k == j)
                BB1[i][j] += L[i][k] * d[k];
            }
        }
    }

    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < 10; j++)
        {
            for (k = 0; k < 10; k++)
            {
                BB2[i][j] += BB1[i][k] * L[j][k];
            }
        }
    }

    printf("BB2\n");
    for (i = 0; i < 10; i++)
    {
        for (j = 0; j < 10; j++)
        {
            printf("%le ", BB2[i][j]);
        }
        printf("\n");
    }
    printf("\n");


    free(A);

}