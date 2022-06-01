#include <iostream>
#include <stdio.h>

using namespace std;

FILE *fp;

int main(int argc, char **argv)
{
    char s[256];
    int temp_i;
    double temp_d;

    fp = fopen(argv[1], "r");

    fscanf(fp, "%lf", &temp_d);
    printf("%le\n", temp_d);

    fgets(s, 256, fp);

    fscanf(fp, "%lf", &temp_d);
    printf("%le\n", temp_d);

    fgets(s, 256, fp);

    fscanf(fp, "%lf", &temp_d);
    printf("%le\n", temp_d);

    fgets(s, 256, fp);

    fscanf(fp, "%d", &temp_i);
    printf("%d\n", temp_i);

    fclose(fp);
}
