/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#include <stdio.h>
#include <stdlib.h>

#include "face.h"

#include "constant.h"
#include "scalar.h"

extern int element_type;
extern int num_elements;
extern int num_element_nodes;
extern int num_element_triangles;
extern int num_element_quadrangles;
extern int **elements;
extern int num_visible_elements;
extern int *visible_elements;
extern int num_triangles;
extern int num_quadrangles;
extern int **triangles;
extern int **quadrangles;
extern int num_visible_triangles;
extern int num_visible_quadrangles;
extern int *visible_triangles;
extern int *visible_quadrangles;
extern int num_dimensions;

static int isSameTrianglePair(int triangle_a_id, int triangle_b_id);
static int isSameQuadranglePair(int quadrangle_a_id, int quadrangle_b_id);
static int isZeroAreaTriangle(int triangle_id);
static int isZeroAreaQuadrangle(int quadrangle_id);

/*
 * Init or reinit faces
 */
void initOrReinitFaces(void)
{
    void *flag;

    static int previous_num_triangles = -1;
    static int previous_num_quadrangles = -1;
    static int previous_num_visible_triangles = -1;
    static int previous_num_visible_quadrangles = -1;

    int i;

    /* Get numbers */
    switch (element_type)
    {
    case LINEAR_TETRAHEDRON_ELEMENT_TYPE:
    case QUADRATIC_TETRAHEDRON_ELEMENT_TYPE:
        num_element_triangles = 4;
        num_element_quadrangles = 0;
        break;
    case LINEAR_HEXAHEDRON_ELEMENT_TYPE:
    case QUADRATIC_HEXAHEDRON_ELEMENT_TYPE:
    case QUADRATIC_HEXAHEDRON27_ELEMENT_TYPE:
        num_element_triangles = 0;
        num_element_quadrangles = 6;
        break;
    case LINEAR_TRIANGLE_ELEMENT_TYPE:
    case QUADRATIC_TRIANGLE_ELEMENT_TYPE:
        num_element_triangles = 1;
        num_element_quadrangles = 0;
        break;
    case LINEAR_QUADRANGLE_ELEMENT_TYPE:
    case QUADRATIC_QUADRANGLE_ELEMENT_TYPE:
    case QUADRATIC_QUADRANGLE9_ELEMENT_TYPE:
        num_element_triangles = 0;
        num_element_quadrangles = 1;
        break;
    default:
        num_element_triangles = 0;
        num_element_quadrangles = 0;
    }
    num_triangles = num_element_triangles * num_elements;
    num_quadrangles = num_element_quadrangles * num_elements;

    /* Allocate faces */
    if (num_triangles > previous_num_triangles)
    {
        flag = (void *)triangles;
        triangles
            = (int **)realloc(triangles,
                              sizeof(int *) * num_triangles);
        if (num_triangles)
        {
            if (flag == NULL) triangles[0] = NULL;
            triangles[0]
                = (int *)realloc(triangles[0],
                                 sizeof(int) * 6 * num_triangles);
        }
        for (i = 0; i < num_triangles; i++)
            triangles[i]
                = triangles[0] + 6 * i;

        previous_num_triangles = num_triangles;
    }
    if (num_quadrangles > previous_num_quadrangles)
    {
        flag = (void *)quadrangles;
        quadrangles
            = (int **)realloc(quadrangles,
                              sizeof(int *) * num_quadrangles);
        if (num_quadrangles)
        {
            if (flag == NULL) quadrangles[0] = NULL;
            quadrangles[0]
                = (int *)realloc(quadrangles[0],
                                 sizeof(int) * 9 * num_quadrangles);
        }
        for (i = 0; i < num_quadrangles; i++)
            quadrangles[i]
                = quadrangles[0] + 9 * i;

        previous_num_quadrangles = num_quadrangles;
    }

    /* Generate faces */
    switch (element_type)
    {
    case LINEAR_TETRAHEDRON_ELEMENT_TYPE:
    case QUADRATIC_TETRAHEDRON_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            triangles[4 * i    ][0] = elements[i][1];
            triangles[4 * i    ][1] = elements[i][2];
            triangles[4 * i    ][2] = elements[i][3];

            triangles[4 * i + 1][0] = elements[i][0];
            triangles[4 * i + 1][1] = elements[i][3];
            triangles[4 * i + 1][2] = elements[i][2];

            triangles[4 * i + 2][0] = elements[i][0];
            triangles[4 * i + 2][1] = elements[i][1];
            triangles[4 * i + 2][2] = elements[i][3];

            triangles[4 * i + 3][0] = elements[i][0];
            triangles[4 * i + 3][1] = elements[i][2];
            triangles[4 * i + 3][2] = elements[i][1];
        }
        break;
    case LINEAR_HEXAHEDRON_ELEMENT_TYPE:
    case QUADRATIC_HEXAHEDRON_ELEMENT_TYPE:
    case QUADRATIC_HEXAHEDRON27_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            quadrangles[6 * i    ][0] = elements[i][0];
            quadrangles[6 * i    ][1] = elements[i][3];
            quadrangles[6 * i    ][2] = elements[i][2];
            quadrangles[6 * i    ][3] = elements[i][1];

            quadrangles[6 * i + 1][0] = elements[i][4];
            quadrangles[6 * i + 1][1] = elements[i][5];
            quadrangles[6 * i + 1][2] = elements[i][6];
            quadrangles[6 * i + 1][3] = elements[i][7];

            quadrangles[6 * i + 2][0] = elements[i][0];
            quadrangles[6 * i + 2][1] = elements[i][1];
            quadrangles[6 * i + 2][2] = elements[i][5];
            quadrangles[6 * i + 2][3] = elements[i][4];

            quadrangles[6 * i + 3][0] = elements[i][1];
            quadrangles[6 * i + 3][1] = elements[i][2];
            quadrangles[6 * i + 3][2] = elements[i][6];
            quadrangles[6 * i + 3][3] = elements[i][5];

            quadrangles[6 * i + 4][0] = elements[i][2];
            quadrangles[6 * i + 4][1] = elements[i][3];
            quadrangles[6 * i + 4][2] = elements[i][7];
            quadrangles[6 * i + 4][3] = elements[i][6];

            quadrangles[6 * i + 5][0] = elements[i][0];
            quadrangles[6 * i + 5][1] = elements[i][4];
            quadrangles[6 * i + 5][2] = elements[i][7];
            quadrangles[6 * i + 5][3] = elements[i][3];
        }
        break;
    case LINEAR_TRIANGLE_ELEMENT_TYPE:
    case QUADRATIC_TRIANGLE_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            triangles[i][0] = elements[i][0];
            triangles[i][1] = elements[i][1];
            triangles[i][2] = elements[i][2];
        }
        break;
    case LINEAR_QUADRANGLE_ELEMENT_TYPE:
    case QUADRATIC_QUADRANGLE_ELEMENT_TYPE:
    case QUADRATIC_QUADRANGLE9_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            quadrangles[i][0] = elements[i][0];
            quadrangles[i][1] = elements[i][1];
            quadrangles[i][2] = elements[i][2];
            quadrangles[i][3] = elements[i][3];
        }
        break;
    }
    switch (element_type)
    {
    case LINEAR_TETRAHEDRON_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            triangles[4 * i    ][3] = -1;
            triangles[4 * i    ][4] = -1;
            triangles[4 * i    ][5] = -1;

            triangles[4 * i + 1][3] = -1;
            triangles[4 * i + 1][4] = -1;
            triangles[4 * i + 1][5] = -1;

            triangles[4 * i + 2][3] = -1;
            triangles[4 * i + 2][4] = -1;
            triangles[4 * i + 2][5] = -1;

            triangles[4 * i + 3][3] = -1;
            triangles[4 * i + 3][4] = -1;
            triangles[4 * i + 3][5] = -1;
        }
        break;
    case QUADRATIC_TETRAHEDRON_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            triangles[4 * i    ][3] = elements[i][7];
            triangles[4 * i    ][4] = elements[i][8];
            triangles[4 * i    ][5] = elements[i][9];

            triangles[4 * i + 1][3] = elements[i][6];
            triangles[4 * i + 1][4] = elements[i][8];
            triangles[4 * i + 1][5] = elements[i][5];

            triangles[4 * i + 2][3] = elements[i][4];
            triangles[4 * i + 2][4] = elements[i][9];
            triangles[4 * i + 2][5] = elements[i][6];

            triangles[4 * i + 3][3] = elements[i][5];
            triangles[4 * i + 3][4] = elements[i][7];
            triangles[4 * i + 3][5] = elements[i][4];
        }
        break;
    case LINEAR_HEXAHEDRON_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            quadrangles[6 * i    ][4] = -1;
            quadrangles[6 * i    ][5] = -1;
            quadrangles[6 * i    ][6] = -1;
            quadrangles[6 * i    ][7] = -1;
            quadrangles[6 * i    ][8] = -1;

            quadrangles[6 * i + 1][4] = -1;
            quadrangles[6 * i + 1][5] = -1;
            quadrangles[6 * i + 1][6] = -1;
            quadrangles[6 * i + 1][7] = -1;
            quadrangles[6 * i + 1][8] = -1;

            quadrangles[6 * i + 2][4] = -1;
            quadrangles[6 * i + 2][5] = -1;
            quadrangles[6 * i + 2][6] = -1;
            quadrangles[6 * i + 2][7] = -1;
            quadrangles[6 * i + 2][8] = -1;

            quadrangles[6 * i + 3][4] = -1;
            quadrangles[6 * i + 3][5] = -1;
            quadrangles[6 * i + 3][6] = -1;
            quadrangles[6 * i + 3][7] = -1;
            quadrangles[6 * i + 3][8] = -1;

            quadrangles[6 * i + 4][4] = -1;
            quadrangles[6 * i + 4][5] = -1;
            quadrangles[6 * i + 4][6] = -1;
            quadrangles[6 * i + 4][7] = -1;
            quadrangles[6 * i + 4][8] = -1;

            quadrangles[6 * i + 5][4] = -1;
            quadrangles[6 * i + 5][5] = -1;
            quadrangles[6 * i + 5][6] = -1;
            quadrangles[6 * i + 5][7] = -1;
            quadrangles[6 * i + 5][8] = -1;
        }
        break;
    case QUADRATIC_HEXAHEDRON_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            quadrangles[6 * i    ][4] = elements[i][11];
            quadrangles[6 * i    ][5] = elements[i][10];
            quadrangles[6 * i    ][6] = elements[i][9];
            quadrangles[6 * i    ][7] = elements[i][8];
            quadrangles[6 * i    ][8] = -1;

            quadrangles[6 * i + 1][4] = elements[i][16];
            quadrangles[6 * i + 1][5] = elements[i][17];
            quadrangles[6 * i + 1][6] = elements[i][18];
            quadrangles[6 * i + 1][7] = elements[i][19];
            quadrangles[6 * i + 1][8] = -1;

            quadrangles[6 * i + 2][4] = elements[i][8];
            quadrangles[6 * i + 2][5] = elements[i][13];
            quadrangles[6 * i + 2][6] = elements[i][16];
            quadrangles[6 * i + 2][7] = elements[i][12];
            quadrangles[6 * i + 2][8] = -1;

            quadrangles[6 * i + 3][4] = elements[i][9];
            quadrangles[6 * i + 3][5] = elements[i][14];
            quadrangles[6 * i + 3][6] = elements[i][17];
            quadrangles[6 * i + 3][7] = elements[i][13];
            quadrangles[6 * i + 3][8] = -1;

            quadrangles[6 * i + 4][4] = elements[i][10];
            quadrangles[6 * i + 4][5] = elements[i][15];
            quadrangles[6 * i + 4][6] = elements[i][18];
            quadrangles[6 * i + 4][7] = elements[i][14];
            quadrangles[6 * i + 4][8] = -1;

            quadrangles[6 * i + 5][4] = elements[i][12];
            quadrangles[6 * i + 5][5] = elements[i][19];
            quadrangles[6 * i + 5][6] = elements[i][15];
            quadrangles[6 * i + 5][7] = elements[i][11];
            quadrangles[6 * i + 5][8] = -1;
        }
        break;
    case QUADRATIC_HEXAHEDRON27_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            quadrangles[6 * i    ][4] = elements[i][11];
            quadrangles[6 * i    ][5] = elements[i][10];
            quadrangles[6 * i    ][6] = elements[i][9];
            quadrangles[6 * i    ][7] = elements[i][8];
            quadrangles[6 * i    ][8] = elements[i][21];

            quadrangles[6 * i + 1][4] = elements[i][16];
            quadrangles[6 * i + 1][5] = elements[i][17];
            quadrangles[6 * i + 1][6] = elements[i][18];
            quadrangles[6 * i + 1][7] = elements[i][19];
            quadrangles[6 * i + 1][8] = elements[i][22];

            quadrangles[6 * i + 2][4] = elements[i][8];
            quadrangles[6 * i + 2][5] = elements[i][13];
            quadrangles[6 * i + 2][6] = elements[i][16];
            quadrangles[6 * i + 2][7] = elements[i][12];
            quadrangles[6 * i + 2][8] = elements[i][23];

            quadrangles[6 * i + 3][4] = elements[i][9];
            quadrangles[6 * i + 3][5] = elements[i][14];
            quadrangles[6 * i + 3][6] = elements[i][17];
            quadrangles[6 * i + 3][7] = elements[i][13];
            quadrangles[6 * i + 3][8] = elements[i][24];

            quadrangles[6 * i + 4][4] = elements[i][10];
            quadrangles[6 * i + 4][5] = elements[i][15];
            quadrangles[6 * i + 4][6] = elements[i][18];
            quadrangles[6 * i + 4][7] = elements[i][14];
            quadrangles[6 * i + 4][8] = elements[i][25];

            quadrangles[6 * i + 5][4] = elements[i][12];
            quadrangles[6 * i + 5][5] = elements[i][19];
            quadrangles[6 * i + 5][6] = elements[i][15];
            quadrangles[6 * i + 5][7] = elements[i][11];
            quadrangles[6 * i + 5][8] = elements[i][26];
        }
        break;
    case LINEAR_TRIANGLE_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            triangles[i][3] = -1;
            triangles[i][4] = -1;
            triangles[i][5] = -1;
        }
        break;
    case QUADRATIC_TRIANGLE_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            triangles[i][3] = elements[i][3];
            triangles[i][4] = elements[i][4];
            triangles[i][5] = elements[i][5];
        }
        break;
    case LINEAR_QUADRANGLE_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            quadrangles[i][4] = -1;
            quadrangles[i][5] = -1;
            quadrangles[i][6] = -1;
            quadrangles[i][7] = -1;
            quadrangles[i][8] = -1;
        }
        break;
    case QUADRATIC_QUADRANGLE_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            quadrangles[i][4] = elements[i][4];
            quadrangles[i][5] = elements[i][5];
            quadrangles[i][6] = elements[i][6];
            quadrangles[i][7] = elements[i][7];
            quadrangles[i][8] = -1;
        }
        break;
    case QUADRATIC_QUADRANGLE9_ELEMENT_TYPE:
        for (i = 0; i < num_elements; i++)
        {
            quadrangles[i][4] = elements[i][4];
            quadrangles[i][5] = elements[i][5];
            quadrangles[i][6] = elements[i][6];
            quadrangles[i][7] = elements[i][7];
            quadrangles[i][8] = elements[i][8];
        }
        break;
    }

    /* Get number of visible faces */
    num_visible_triangles   = num_triangles;
    num_visible_quadrangles = num_quadrangles;

    /* Allocate visible faces */
    if (num_visible_triangles > previous_num_visible_triangles)
    {
        visible_triangles
            = (int *)realloc(visible_triangles,
                             sizeof(int) * num_visible_triangles);

        previous_num_visible_triangles = num_visible_triangles;
    }
    if (num_visible_quadrangles > previous_num_visible_quadrangles)
    {
        visible_quadrangles
            = (int *)realloc(visible_quadrangles,
                             sizeof(int) * num_visible_quadrangles);

        previous_num_visible_quadrangles = num_visible_quadrangles;
    }

    /* Fill visible faces */
    for (i = 0; i < num_visible_triangles; i++)
        visible_triangles[i] = i;
    for (i = 0; i < num_visible_quadrangles; i++)
        visible_quadrangles[i] = i;

    /* Remove inner faces */
    removeInnerTriangles();
    removeInnerQuadrangles();
}

/*
 * Destroy faces
 */
void destroyFaces(void)
{
    free(visible_quadrangles);
    free(visible_triangles);

    if (num_quadrangles)
        free(quadrangles[0]);
    free(quadrangles);

    if (num_triangles)
        free(triangles[0]);
    free(triangles);
}

/*
 * Reinit faces from visible elements
 */
void reinitFacesFromVisibleElements(void)
{
    int i, j;

    /* Get number of visible faces */
    num_visible_triangles = num_element_triangles * num_visible_elements;
    num_visible_quadrangles = num_element_quadrangles * num_visible_elements;

    /* Generate visible faces */
    for (i = 0; i < num_visible_elements; i++)
        for (j = 0; j < num_element_triangles; j++)
            visible_triangles[num_element_triangles * i + j]
                = num_element_triangles * visible_elements[i] + j;
    for (i = 0; i < num_visible_elements; i++)
        for (j = 0; j < num_element_quadrangles; j++)
            visible_quadrangles[num_element_quadrangles * i + j]
                = num_element_quadrangles * visible_elements[i] + j;

    /* Remove inner faces */
    removeInnerTriangles();
    removeInnerQuadrangles();
}

/*
 * Remove inner triangles
 */
void removeInnerTriangles(void)
{
    int *hashes;
    int max_hash;
    double hash_factor;

    int *hash_table_ind;
    int *hash_table_ptr;
    int hash_table_size = num_visible_triangles;

    int offset;
    int i, j, k;

    if (num_dimensions != 3)
        return;

    if (num_visible_triangles == 0)
        return;

    /* Allocate hashes and hash table */
    hashes = (int *)malloc(sizeof(int) * num_visible_triangles);
    hash_table_ind = (int *)malloc(sizeof(int) * num_visible_triangles);
    hash_table_ptr = (int *)malloc(sizeof(int) * (hash_table_size + 1));

    /* Calculate hashes */
    max_hash = 1;
    for (i = 0; i < num_visible_triangles; i++)
    {
        int triangle_id = visible_triangles[i];

        int node0_id = triangles[triangle_id][0];
        int node1_id = triangles[triangle_id][1];
        int node2_id = triangles[triangle_id][2];

        sortThreeIntegers(&node0_id, &node1_id, &node2_id);

        if (node0_id == node1_id) node0_id = 0;
        if (node1_id == node2_id) node1_id = 0;

        hashes[i]
            = node0_id
            + node1_id
            + node2_id;

        if (max_hash < hashes[i])
            max_hash = hashes[i];
    }
    hash_factor = (double)(hash_table_size - 1) / (double)max_hash;
    for (i = 0; i < num_visible_triangles; i++)
        hashes[i] = (int)(hash_factor * (double)hashes[i]);

    /* Generate hash table */
    for (i = 0; i < hash_table_size + 1; i++)
        hash_table_ptr[i] = 0;
    for (i = 0; i < num_visible_triangles; i++)
        hash_table_ptr[hashes[i] + 1]++;
    for (i = 0; i < hash_table_size; i++)
        hash_table_ptr[i + 1] += hash_table_ptr[i];
    for (i = 0; i < num_visible_triangles; i++)
        hash_table_ind[i] = -1;
    for (i = 0; i < num_visible_triangles; i++)
        for (j = hash_table_ptr[hashes[i]];
             j < hash_table_ptr[hashes[i] + 1];
             j++)
            if (hash_table_ind[j] == -1)
            {
                hash_table_ind[j] = i;
                break;
            }

    /* Find visible triangles */
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
            if (isZeroAreaTriangle(visible_triangles[hash_table_ind[j]]))
                visible_triangles[hash_table_ind[j]] = -1;
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (visible_triangles[hash_table_ind[j]] == -1) continue;

            for (k = j + 1; k < hash_table_ptr[i + 1]; k++)
            {
                if (visible_triangles[hash_table_ind[k]] == -1) continue;

                if (isSameTrianglePair(visible_triangles[hash_table_ind[j]],
                                       visible_triangles[hash_table_ind[k]]))
                {
                    visible_triangles[hash_table_ind[j]] = -1;
                    visible_triangles[hash_table_ind[k]] = -1;
                    break;
                }
            }
        }
    offset = 0;
    for (i = 0; i < num_visible_triangles; i++)
        if (visible_triangles[i] != -1)
            visible_triangles[offset++] = visible_triangles[i];
    num_visible_triangles = offset;

    /* Destroy hahses and hash table */
    free(hash_table_ptr);
    free(hash_table_ind);
    free(hashes);
}

/*
 * remove inner quadrangles
 */
void removeInnerQuadrangles(void)
{
    int *hashes;
    int max_hash;
    double hash_factor;

    int *hash_table_ind;
    int *hash_table_ptr;
    int hash_table_size = num_visible_quadrangles;

    int offset;
    int i, j, k;

    if (num_dimensions != 3)
        return;

    if (num_visible_quadrangles == 0)
        return;

    /* Allocate hashes and hash table */
    hashes = (int *)malloc(sizeof(int) * num_visible_quadrangles);
    hash_table_ind = (int *)malloc(sizeof(int) * num_visible_quadrangles);
    hash_table_ptr = (int *)malloc(sizeof(int) * (hash_table_size + 1));

    /* Calculate hashes */
    max_hash = 1;
    for (i = 0; i < num_visible_quadrangles; i++)
    {
        int quadrangle_id = visible_quadrangles[i];

        int node0_id = quadrangles[quadrangle_id][0];
        int node1_id = quadrangles[quadrangle_id][1];
        int node2_id = quadrangles[quadrangle_id][2];
        int node3_id = quadrangles[quadrangle_id][3];

        sortFourIntegers(&node0_id, &node1_id, &node2_id, &node3_id);

        if (node0_id == node1_id) node0_id = 0;
        if (node1_id == node2_id) node1_id = 0;
        if (node2_id == node3_id) node2_id = 0;

        hashes[i]
            = node0_id
            + node1_id
            + node2_id
            + node3_id;

        if (max_hash < hashes[i])
            max_hash = hashes[i];
    }
    hash_factor = (double)(hash_table_size - 1) / (double)max_hash;
    for (i = 0; i < num_visible_quadrangles; i++)
        hashes[i] = (int)(hash_factor * (double)hashes[i]);

    /* Generate hash table */
    for (i = 0; i < hash_table_size + 1; i++)
        hash_table_ptr[i] = 0;
    for (i = 0; i < num_visible_quadrangles; i++)
        hash_table_ptr[hashes[i] + 1]++;
    for (i = 0; i < hash_table_size; i++)
        hash_table_ptr[i + 1] += hash_table_ptr[i];
    for (i = 0; i < num_visible_quadrangles; i++)
        hash_table_ind[i] = -1;
    for (i = 0; i < num_visible_quadrangles; i++)
        for (j = hash_table_ptr[hashes[i]];
             j < hash_table_ptr[hashes[i] + 1];
             j++)
            if (hash_table_ind[j] == -1)
            {
                hash_table_ind[j] = i;
                break;
            }

    /* Find visible quadrangles */
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
            if (isZeroAreaQuadrangle(visible_quadrangles[hash_table_ind[j]]))
                visible_quadrangles[hash_table_ind[j]] = -1;
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (visible_quadrangles[hash_table_ind[j]] == -1) continue;

            for (k = j + 1; k < hash_table_ptr[i + 1]; k++)
            {
                if (visible_quadrangles[hash_table_ind[k]] == -1) continue;

                if (isSameQuadranglePair(visible_quadrangles[hash_table_ind[j]],
                                         visible_quadrangles[hash_table_ind[k]]))
                {
                    visible_quadrangles[hash_table_ind[j]] = -1;
                    visible_quadrangles[hash_table_ind[k]] = -1;
                    break;
                }
            }
        }
    offset = 0;
    for (i = 0; i < num_visible_quadrangles; i++)
        if (visible_quadrangles[i] != -1)
            visible_quadrangles[offset++] = visible_quadrangles[i];
    num_visible_quadrangles = offset;

    /* Destroy hahses and hash table */
    free(hash_table_ptr);
    free(hash_table_ind);
    free(hashes);
}

/*
 * Is same triangle pair?
 */
static int isSameTrianglePair(int triangle_a_id, int triangle_b_id)
{
    int min_a = triangles[triangle_a_id][0];
    int mid_a = triangles[triangle_a_id][1];
    int max_a = triangles[triangle_a_id][2];
    int min_b = triangles[triangle_b_id][0];
    int mid_b = triangles[triangle_b_id][1];
    int max_b = triangles[triangle_b_id][2];

    sortThreeIntegers(&min_a, &mid_a, &max_a);
    sortThreeIntegers(&min_b, &mid_b, &max_b);

    if (min_a == mid_a) min_a = -1;
    if (mid_a == max_a) mid_a = -1;
    if (min_b == mid_b) min_b = -1;
    if (mid_b == max_b) mid_b = -1;

    sortThreeIntegers(&min_a, &mid_a, &max_a);
    sortThreeIntegers(&min_b, &mid_b, &max_b);

    return
        min_a == min_b
        && mid_a == mid_b
        && max_a == max_b;
}

/*
 * Is same quadrangle pair?
 */
static int isSameQuadranglePair(int quadrangle_a_id, int quadrangle_b_id)
{
    int minmin_a = quadrangles[quadrangle_a_id][0];
    int midmin_a = quadrangles[quadrangle_a_id][1];
    int midmax_a = quadrangles[quadrangle_a_id][2];
    int maxmax_a = quadrangles[quadrangle_a_id][3];
    int minmin_b = quadrangles[quadrangle_b_id][0];
    int midmin_b = quadrangles[quadrangle_b_id][1];
    int midmax_b = quadrangles[quadrangle_b_id][2];
    int maxmax_b = quadrangles[quadrangle_b_id][3];

    sortFourIntegers(&minmin_a, &midmin_a, &midmax_a, &maxmax_a);
    sortFourIntegers(&minmin_b, &midmin_b, &midmax_b, &maxmax_b);

    if (minmin_a == midmin_a) minmin_a = -1;
    if (midmin_a == midmax_a) midmin_a = -1;
    if (midmax_a == maxmax_a) midmax_a = -1;
    if (minmin_b == midmin_b) minmin_b = -1;
    if (midmin_b == midmax_b) midmin_b = -1;
    if (midmax_b == maxmax_b) midmax_b = -1;

    sortFourIntegers(&minmin_a, &midmin_a, &midmax_a, &maxmax_a);
    sortFourIntegers(&minmin_b, &midmin_b, &midmax_b, &maxmax_b);

    return
        minmin_a == minmin_b
        && midmin_a == midmin_b
        && midmax_a == midmax_b
        && maxmax_a == maxmax_b;
}

/*
 * Is zero-area triangle?
 */
static int isZeroAreaTriangle(int triangle_id)
{
    int min = triangles[triangle_id][0];
    int mid = triangles[triangle_id][1];
    int max = triangles[triangle_id][2];

    sortThreeIntegers(&min, &mid, &max);

    return
        min == mid
        || mid == max;
}

/*
 * Is zero-area quadrangle?
 */
static int isZeroAreaQuadrangle(int quadrangle_id)
{
    int minmin = quadrangles[quadrangle_id][0];
    int midmin = quadrangles[quadrangle_id][1];
    int midmax = quadrangles[quadrangle_id][2];
    int maxmax = quadrangles[quadrangle_id][3];

    sortFourIntegers(&minmin, &midmin, &midmax, &maxmax);

    return
        (minmin == midmin && midmax == maxmax)
        || (minmin == midmin && midmin == midmax)
        || (midmin == midmax && midmax == maxmax);
}
