/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "crack.h"

#include "constant.h"
#include "edge.h"
#include "scalar.h"

extern int **triangles;
extern int **quadrangles;
extern int num_visible_triangles;
extern int num_visible_quadrangles;
extern int *visible_triangles;
extern int *visible_quadrangles;
extern int num_crack_triangles;
extern int num_crack_quadrangles;
extern int *crack_triangles;
extern int *crack_quadrangles;
extern int num_crack_front_triangle_edges;
extern int num_crack_front_quadrangle_edges;
extern int *crack_front_triangle_edges;
extern int *crack_front_quadrangle_edges;
extern int num_dimensions;
extern double **nodes;

static int isCrackTrianglePair(int triangle_a_id, int triangle_b_id);
static int isCrackQuadranglePair(int quadrangle_a_id, int quadrangle_b_id);
static int compareNodalCoords(const void *a, const void *b);

/*
 * Init or reinit crack faces
 */
void initOrReinitCrackFaces(void)
{
    static int previous_num_crack_triangles = -1;
    static int previous_num_crack_quadrangles = -1;

    int i;

    if (num_dimensions != 3)
        return;

    /* Get number of crack faces */
    num_crack_triangles   = num_visible_triangles;
    num_crack_quadrangles = num_visible_quadrangles;

    /* Allocate crack faces */
    if (num_crack_triangles > previous_num_crack_triangles)
    {
        crack_triangles
            = (int *)realloc(crack_triangles,
                             sizeof(int) * num_crack_triangles);

        previous_num_crack_triangles = num_crack_triangles;
    }
    if (num_crack_quadrangles > previous_num_crack_quadrangles)
    {
        crack_quadrangles
            = (int *)realloc(crack_quadrangles,
                             sizeof(int) * num_crack_quadrangles);

        previous_num_crack_quadrangles = num_crack_quadrangles;
    }

    /* Fill crack faces */
    for (i = 0; i < num_crack_triangles; i++)
        crack_triangles[i] = visible_triangles[i];
    for (i = 0; i < num_crack_quadrangles; i++)
        crack_quadrangles[i] = visible_quadrangles[i];

    /* Remove non-crack faces */
    removeNonCrackTriangles();
    removeNonCrackQuadrangles();
}

/*
 * Destroy crack faces
 */
void destroyCrackFaces(void)
{
    if (num_dimensions != 3)
        return;

    free(crack_quadrangles);
    free(crack_triangles);
}

/*
 * Remove non-crack triangles
 */
void removeNonCrackTriangles(void)
{
    int *hashes;
    double *hashes2;
    double max_hash, min_hash;
    double hash_factor;

    int *hash_table_ind;
    int *hash_table_ptr;
    int hash_table_size = num_crack_triangles;

    int offset;
    int i, j, k;

    if (num_crack_triangles == 0) return;

    /* Allocate hashes and hash table */
    hashes = (int *)malloc(sizeof(int) * num_crack_triangles);
    hashes2 = (double *)malloc(sizeof(double) * num_crack_triangles);
    hash_table_ind = (int *)malloc(sizeof(int) * num_crack_triangles);
    hash_table_ptr = (int *)malloc(sizeof(int) * (hash_table_size + 1));

    /* Calculate hashes */
    max_hash = -DBL_MAX;
    min_hash = DBL_MAX;
    for (i = 0; i < num_crack_triangles; i++)
    {
        int triangle_id = crack_triangles[i];

        int node0_id = triangles[triangle_id][0];
        int node1_id = triangles[triangle_id][1];
        int node2_id = triangles[triangle_id][2];

        double x0, y0, z0;
        double x1, y1, z1;
        double x2, y2, z2;

        sortThreeIntegers(&node0_id, &node1_id, &node2_id);

        x0 = nodes[node0_id][0];
        y0 = nodes[node0_id][1];
        z0 = nodes[node0_id][2];
        x1 = nodes[node1_id][0];
        y1 = nodes[node1_id][1];
        z1 = nodes[node1_id][2];
        x2 = nodes[node2_id][0];
        y2 = nodes[node2_id][1];
        z2 = nodes[node2_id][2];

        if (node0_id == node1_id) {x0 = 0.0; y0 = 0.0; z0 = 0.0;}
        if (node1_id == node2_id) {x1 = 0.0; y1 = 0.0; z1 = 0.0;}

        hashes2[i]
            = x0 + y0 + z0
            + x1 + y1 + z1
            + x2 + y2 + z2;

        if (max_hash < hashes2[i])
            max_hash = hashes2[i];
        if (min_hash > hashes2[i])
            min_hash = hashes2[i];
    }
    if (max_hash > min_hash)
        hash_factor = (double)(hash_table_size - 1) / (max_hash - min_hash);
    else
        hash_factor = 0.0;
    for (i = 0; i < num_crack_triangles; i++)
        hashes[i] = (int)(hash_factor * (hashes2[i] - min_hash));

    /* Generate hash table */
    for (i = 0; i < hash_table_size + 1; i++)
        hash_table_ptr[i] = 0;
    for (i = 0; i < num_crack_triangles; i++)
        hash_table_ptr[hashes[i] + 1]++;
    for (i = 0; i < hash_table_size; i++)
        hash_table_ptr[i + 1] += hash_table_ptr[i];
    for (i = 0; i < num_crack_triangles; i++)
        hash_table_ind[i] = -1;
    for (i = 0; i < num_crack_triangles; i++)
        for (j = hash_table_ptr[hashes[i]];
             j < hash_table_ptr[hashes[i] + 1];
             j++)
            if (hash_table_ind[j] == -1)
            {
                hash_table_ind[j] = i;
                break;
            }

    /* Find crack triangles */
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (crack_triangles[hash_table_ind[j]] < 0) continue;

            for (k = j + 1; k < hash_table_ptr[(i < hash_table_size - 1)? i + 2: i + 1]; k++)
            {
                if (crack_triangles[hash_table_ind[k]] < 0) continue;

                if (isCrackTrianglePair(crack_triangles[hash_table_ind[j]],
                                        crack_triangles[hash_table_ind[k]]))
                {
                    crack_triangles[hash_table_ind[j]]
                        = -crack_triangles[hash_table_ind[j]] - 1;
                    crack_triangles[hash_table_ind[k]]
                        = -crack_triangles[hash_table_ind[k]] - 1;
                    break;
                }
            }
        }
    offset = 0;
    for (i = 0; i < num_crack_triangles; i++)
        if (crack_triangles[i] < 0)
            crack_triangles[offset++] = -crack_triangles[i] - 1;
    num_crack_triangles = offset;

    /* Destroy hashes and hash table */
    free(hash_table_ptr);
    free(hash_table_ind);
    free(hashes2);
    free(hashes);
}

/*
 * Remove non-crack quadrangles
 */
void removeNonCrackQuadrangles(void)
{
    int *hashes;
    double *hashes2;
    double max_hash, min_hash;
    double hash_factor;

    int *hash_table_ind;
    int *hash_table_ptr;
    int hash_table_size = num_crack_quadrangles;

    int offset;
    int i, j, k;

    if (num_crack_quadrangles == 0) return;

    /* Allocate hashes and hash table */
    hashes = (int *)malloc(sizeof(int) * num_crack_quadrangles);
    hashes2 = (double *)malloc(sizeof(double) * num_crack_quadrangles);
    hash_table_ind = (int *)malloc(sizeof(int) * num_crack_quadrangles);
    hash_table_ptr = (int *)malloc(sizeof(int) * (hash_table_size + 1));

    /* Calculate hashes */
    max_hash = -DBL_MAX;
    min_hash = DBL_MAX;
    for (i = 0; i < num_crack_quadrangles; i++)
    {
        int quadrangle_id = crack_quadrangles[i];

        int node0_id = quadrangles[quadrangle_id][0];
        int node1_id = quadrangles[quadrangle_id][1];
        int node2_id = quadrangles[quadrangle_id][2];
        int node3_id = quadrangles[quadrangle_id][3];

        double x0, y0, z0;
        double x1, y1, z1;
        double x2, y2, z2;
        double x3, y3, z3;

        sortFourIntegers(&node0_id, &node1_id, &node2_id, &node3_id);

        x0 = nodes[node0_id][0];
        y0 = nodes[node0_id][1];
        z0 = nodes[node0_id][2];
        x1 = nodes[node1_id][0];
        y1 = nodes[node1_id][1];
        z1 = nodes[node1_id][2];
        x2 = nodes[node2_id][0];
        y2 = nodes[node2_id][1];
        z2 = nodes[node2_id][2];
        x3 = nodes[node3_id][0];
        y3 = nodes[node3_id][1];
        z3 = nodes[node3_id][2];

        if (node0_id == node1_id) {x0 = 0.0; y0 = 0.0; z0 = 0.0;}
        if (node1_id == node2_id) {x1 = 0.0; y1 = 0.0; z1 = 0.0;}
        if (node2_id == node3_id) {x2 = 0.0; y2 = 0.0; z2 = 0.0;}

        hashes2[i]
            = x0 + y0 + z0
            + x1 + y1 + z1
            + x2 + y2 + z2
            + x3 + y3 + z3;

        if (max_hash < hashes2[i])
            max_hash = hashes2[i];
        if (min_hash > hashes2[i])
            min_hash = hashes2[i];
    }
    if (max_hash > min_hash)
        hash_factor = (double)(hash_table_size - 1) / (max_hash - min_hash);
    else
        hash_factor = 0.0;
    for (i = 0; i < num_crack_quadrangles; i++)
        hashes[i] = (int)(hash_factor * (hashes2[i] - min_hash));

    /* Generate hash table */
    for (i = 0; i < hash_table_size + 1; i++)
        hash_table_ptr[i] = 0;
    for (i = 0; i < num_crack_quadrangles; i++)
        hash_table_ptr[hashes[i] + 1]++;
    for (i = 0; i < hash_table_size; i++)
        hash_table_ptr[i + 1] += hash_table_ptr[i];
    for (i = 0; i < num_crack_quadrangles; i++)
        hash_table_ind[i] = -1;
    for (i = 0; i < num_crack_quadrangles; i++)
        for (j = hash_table_ptr[hashes[i]];
             j < hash_table_ptr[hashes[i] + 1];
             j++)
            if (hash_table_ind[j] == -1)
            {
                hash_table_ind[j] = i;
                break;
            }

    /* Find crack quadrangles */
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (crack_quadrangles[hash_table_ind[j]] < 0) continue;

            for (k = j + 1; k < hash_table_ptr[(i < hash_table_size - 1)? i + 2: i + 1]; k++)
            {
                if (crack_quadrangles[hash_table_ind[k]] < 0) continue;

                if (isCrackQuadranglePair(crack_quadrangles[hash_table_ind[j]],
                                          crack_quadrangles[hash_table_ind[k]]))
                {
                    crack_quadrangles[hash_table_ind[j]]
                        = -crack_quadrangles[hash_table_ind[j]] - 1;
                    crack_quadrangles[hash_table_ind[k]]
                        = -crack_quadrangles[hash_table_ind[k]] - 1;
                    break;
                }
            }
        }
    offset = 0;
    for (i = 0; i < num_crack_quadrangles; i++)
        if (crack_quadrangles[i] < 0)
            crack_quadrangles[offset++] = -crack_quadrangles[i] - 1;
    num_crack_quadrangles = offset;

    /* Destroy hashes and hash table */
    free(hash_table_ptr);
    free(hash_table_ind);
    free(hashes2);
    free(hashes);
}

/*
 * Is crack triangle pair?
 */
static int isCrackTrianglePair(int triangle_a_id, int triangle_b_id)
{
    int min_a = triangles[triangle_a_id][0];
    int mid_a = triangles[triangle_a_id][1];
    int max_a = triangles[triangle_a_id][2];
    int min_b = triangles[triangle_b_id][0];
    int mid_b = triangles[triangle_b_id][1];
    int max_b = triangles[triangle_b_id][2];

    sortThreeIntegersWithFunction(&min_a, &mid_a, &max_a,
                                  compareNodalCoords);
    sortThreeIntegersWithFunction(&min_b, &mid_b, &max_b,
                                  compareNodalCoords);

    if (min_a == mid_a) min_a = max_a;
    if (mid_a == max_a) mid_a = max_a;
    if (min_b == mid_b) min_b = max_b;
    if (mid_b == max_b) mid_b = max_b;

    sortThreeIntegersWithFunction(&min_a, &mid_a, &max_a,
                                  compareNodalCoords);
    sortThreeIntegersWithFunction(&min_b, &mid_b, &max_b,
                                  compareNodalCoords);

    return
        compareNodalCoords(&min_a, &min_b) == 0
        && compareNodalCoords(&mid_a, &mid_b) == 0
        && compareNodalCoords(&max_a, &max_b) == 0;
}

/*
 * Is crack quadrangle pair?
 */
static int isCrackQuadranglePair(int quadrangle_a_id, int quadrangle_b_id)
{
    int minmin_a = quadrangles[quadrangle_a_id][0];
    int midmin_a = quadrangles[quadrangle_a_id][1];
    int midmax_a = quadrangles[quadrangle_a_id][2];
    int maxmax_a = quadrangles[quadrangle_a_id][3];
    int minmin_b = quadrangles[quadrangle_b_id][0];
    int midmin_b = quadrangles[quadrangle_b_id][1];
    int midmax_b = quadrangles[quadrangle_b_id][2];
    int maxmax_b = quadrangles[quadrangle_b_id][3];

    sortFourIntegersWithFunction(&minmin_a, &midmin_a, &midmax_a, &maxmax_a,
                                 compareNodalCoords);
    sortFourIntegersWithFunction(&minmin_b, &midmin_b, &midmax_b, &maxmax_b,
                                 compareNodalCoords);

    if (minmin_a == midmin_a) minmin_a = maxmax_a;
    if (midmin_a == midmax_a) midmin_a = maxmax_a;
    if (midmax_a == maxmax_a) midmax_a = maxmax_a;
    if (minmin_b == midmin_b) minmin_b = maxmax_b;
    if (midmin_b == midmax_b) midmin_b = maxmax_b;
    if (midmax_b == maxmax_b) midmax_b = maxmax_b;

    sortFourIntegersWithFunction(&minmin_a, &midmin_a, &midmax_a, &maxmax_a,
                                 compareNodalCoords);
    sortFourIntegersWithFunction(&minmin_b, &midmin_b, &midmax_b, &maxmax_b,
                                 compareNodalCoords);

    return
        compareNodalCoords(&minmin_a, &minmin_b) == 0
        && compareNodalCoords(&midmin_a, &midmin_b) == 0
        && compareNodalCoords(&midmax_a, &midmax_b) == 0
        && compareNodalCoords(&maxmax_a, &maxmax_b) == 0;
}

/*
 * Compare nodal coordinates
 */
static int compareNodalCoords(const void *a, const void *b)
{
    const int aa = *(const int *)a;
    const int bb = *(const int *)b;

    double x = nodes[aa][0] - nodes[bb][0];
    double y = nodes[aa][1] - nodes[bb][1];
    double z = nodes[aa][2] - nodes[bb][2];

    int x_eq_zero = (fabs(x) <= EFFECTIVE_ZERO);
    int y_eq_zero = (fabs(y) <= EFFECTIVE_ZERO);
    int z_eq_zero = (fabs(z) <= EFFECTIVE_ZERO);

    if (x_eq_zero && y_eq_zero && z_eq_zero)
        return 0;
    else if (x_eq_zero && y_eq_zero)
        return (z > 0.0)? 1: -1;
    else if (x_eq_zero)
        return (y > 0.0)? 1: -1;
    else
        return (x > 0.0)? 1: -1;
}

/*
 * Init or reinit crack front edges
 */
void initOrReinitCrackFrontEdges(void)
{
    static int previous_num_crack_front_triangle_edges = -1;
    static int previous_num_crack_front_quadrangle_edges = -1;

    int i;

    /* Get number of crack front edges */
    num_crack_front_triangle_edges   = 3 * num_crack_triangles;
    num_crack_front_quadrangle_edges = 4 * num_crack_quadrangles;

    /* Allocate crack front edges */
    if (num_crack_front_triangle_edges > previous_num_crack_front_triangle_edges)
    {
        crack_front_triangle_edges
            = (int *)realloc(crack_front_triangle_edges,
                             sizeof(int) * num_crack_front_triangle_edges);

        previous_num_crack_front_triangle_edges = num_crack_front_triangle_edges;
    }
    if (num_crack_front_quadrangle_edges > previous_num_crack_front_quadrangle_edges)
    {
        crack_front_quadrangle_edges
            = (int *)realloc(crack_front_quadrangle_edges,
                             sizeof(int) * num_crack_front_quadrangle_edges);

        previous_num_crack_front_quadrangle_edges = num_crack_front_quadrangle_edges;
    }

    /* Generate crack front edges */
    for (i = 0; i < num_crack_triangles; i++)
    {
        crack_front_triangle_edges[3 * i    ] = 3 * crack_triangles[i];
        crack_front_triangle_edges[3 * i + 1] = 3 * crack_triangles[i] + 1;
        crack_front_triangle_edges[3 * i + 2] = 3 * crack_triangles[i] + 2;
    }
    for (i = 0; i < num_crack_quadrangles; i++)
    {
        crack_front_quadrangle_edges[4 * i    ] = 4 * crack_quadrangles[i];
        crack_front_quadrangle_edges[4 * i + 1] = 4 * crack_quadrangles[i] + 1;
        crack_front_quadrangle_edges[4 * i + 2] = 4 * crack_quadrangles[i] + 2;
        crack_front_quadrangle_edges[4 * i + 3] = 4 * crack_quadrangles[i] + 3;
    }

    /* Remove non-crack-front edges */
    removeNonCrackFrontTriangleEdges();
    removeNonCrackFrontQuadrangleEdges();
}

/*
 * Destroy crack front edges
 */
void destroyCrackFrontEdges(void)
{
    free(crack_front_quadrangle_edges);
    free(crack_front_triangle_edges);
}

/*
 * Remove non-crack-front triangle edges
 */
void removeNonCrackFrontTriangleEdges(void)
{
    int *hashes;
    int max_hash;
    double hash_factor;

    int *hash_table_ind;
    int *hash_table_ptr;
    int hash_table_size = num_crack_front_triangle_edges;

    int offset;
    int i, j, k;

    if (num_crack_front_triangle_edges == 0) return;

    /* Allocate hashes and hash table */
    hashes = (int *)malloc(sizeof(int) * num_crack_front_triangle_edges);
    hash_table_ind = (int *)malloc(sizeof(int) * num_crack_front_triangle_edges);
    hash_table_ptr = (int *)malloc(sizeof(int) * (hash_table_size + 1));

    /* Calculate hashes */
    max_hash = 1;
    for (i = 0; i < num_crack_front_triangle_edges; i++)
    {
        int triangle_id      = crack_front_triangle_edges[i] / 3;
        int triangle_edge_id = crack_front_triangle_edges[i] % 3;

        int node0_id = triangles[triangle_id][triangle_edge_id];
        int node1_id = triangles[triangle_id][(triangle_edge_id + 1) % 3];

        hashes[i]
            = node0_id
            + node1_id;

        if (max_hash < hashes[i])
            max_hash = hashes[i];
    }
    hash_factor = (double)(hash_table_size - 1) / (double)max_hash;
    for (i = 0; i < num_crack_front_triangle_edges; i++)
        hashes[i] = (int)(hash_factor * (double)hashes[i]);

    /* Generate hash table */
    for (i = 0; i < hash_table_size + 1; i++)
        hash_table_ptr[i] = 0;
    for (i = 0; i < num_crack_front_triangle_edges; i++)
        hash_table_ptr[hashes[i] + 1]++;
    for (i = 0; i < hash_table_size; i++)
        hash_table_ptr[i + 1] += hash_table_ptr[i];
    for (i = 0; i < num_crack_front_triangle_edges; i++)
        hash_table_ind[i] = -1;
    for (i = 0; i < num_crack_front_triangle_edges; i++)
        for (j = hash_table_ptr[hashes[i]];
             j < hash_table_ptr[hashes[i] + 1];
             j++)
            if (hash_table_ind[j] == -1)
            {
                hash_table_ind[j] = i;
                break;
            }

    /* Find crack front triangle edges */
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (crack_front_triangle_edges[hash_table_ind[j]] < 0) continue;

            for (k = j + 1; k < hash_table_ptr[i + 1]; k++)
            {
                if (crack_front_triangle_edges[hash_table_ind[k]] < 0) continue;

                if (isSameTriangleEdgePair(crack_front_triangle_edges[hash_table_ind[j]],
                                           crack_front_triangle_edges[hash_table_ind[k]])
                    && isNonSmoothTriangleEdgePair(crack_front_triangle_edges[hash_table_ind[j]],
                                                   crack_front_triangle_edges[hash_table_ind[k]],
                                                   -1.0 + EFFECTIVE_ZERO))
                {
                    crack_front_triangle_edges[hash_table_ind[j]]
                        = -crack_front_triangle_edges[hash_table_ind[j]] - 1;
                    crack_front_triangle_edges[hash_table_ind[k]]
                        = -crack_front_triangle_edges[hash_table_ind[k]] - 1;
                    break;
                }
            }
        }
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (crack_front_triangle_edges[hash_table_ind[j]] >= 0) continue;

            if (isZeroLengthTriangleEdge(-crack_front_triangle_edges[hash_table_ind[j]] - 1))
                crack_front_triangle_edges[hash_table_ind[j]]
                    = -crack_front_triangle_edges[hash_table_ind[j]] - 1;
        }
    offset = 0;
    for (i = 0; i < num_crack_front_triangle_edges; i++)
        if (crack_front_triangle_edges[i] < 0)
            crack_front_triangle_edges[offset++]
                = -crack_front_triangle_edges[i] - 1;
    num_crack_front_triangle_edges = offset;

    /* Destroy hahses and hash table */
    free(hash_table_ptr);
    free(hash_table_ind);
    free(hashes);
}

/*
 * Remove non-crack-front quadrangle edges
 */
void removeNonCrackFrontQuadrangleEdges(void)
{
    int *hashes;
    int max_hash;
    double hash_factor;

    int *hash_table_ind;
    int *hash_table_ptr;
    int hash_table_size = num_crack_front_quadrangle_edges;

    int offset;
    int i, j, k;

    if (num_crack_front_quadrangle_edges == 0) return;

    /* Allocate hashes and hash table */
    hashes = (int *)malloc(sizeof(int) * num_crack_front_quadrangle_edges);
    hash_table_ind = (int *)malloc(sizeof(int) * num_crack_front_quadrangle_edges);
    hash_table_ptr = (int *)malloc(sizeof(int) * (hash_table_size + 1));

    /* Calculate hashes */
    max_hash = 1;
    for (i = 0; i < num_crack_front_quadrangle_edges; i++)
    {
        int quadrangle_id      = crack_front_quadrangle_edges[i] / 4;
        int quadrangle_edge_id = crack_front_quadrangle_edges[i] % 4;

        int node0_id = quadrangles[quadrangle_id][quadrangle_edge_id];
        int node1_id = quadrangles[quadrangle_id][(quadrangle_edge_id + 1) % 4];

        hashes[i]
            = node0_id
            + node1_id;

        if (max_hash < hashes[i])
            max_hash = hashes[i];
    }
    hash_factor = (double)(hash_table_size - 1) / (double)max_hash;
    for (i = 0; i < num_crack_front_quadrangle_edges; i++)
        hashes[i] = (int)(hash_factor * (double)hashes[i]);

    /* Generate hash table */
    for (i = 0; i < hash_table_size + 1; i++)
        hash_table_ptr[i] = 0;
    for (i = 0; i < num_crack_front_quadrangle_edges; i++)
        hash_table_ptr[hashes[i] + 1]++;
    for (i = 0; i < hash_table_size; i++)
        hash_table_ptr[i + 1] += hash_table_ptr[i];
    for (i = 0; i < num_crack_front_quadrangle_edges; i++)
        hash_table_ind[i] = -1;
    for (i = 0; i < num_crack_front_quadrangle_edges; i++)
        for (j = hash_table_ptr[hashes[i]];
             j < hash_table_ptr[hashes[i] + 1];
             j++)
            if (hash_table_ind[j] == -1)
            {
                hash_table_ind[j] = i;
                break;
            }

    /* Find crack front quadrangle edges */
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (crack_front_quadrangle_edges[hash_table_ind[j]] < 0) continue;

            for (k = j + 1; k < hash_table_ptr[i + 1]; k++)
            {
                if (crack_front_quadrangle_edges[hash_table_ind[k]] < 0) continue;

                if (isSameQuadrangleEdgePair(crack_front_quadrangle_edges[hash_table_ind[j]],
                                           crack_front_quadrangle_edges[hash_table_ind[k]])
                    && isNonSmoothQuadrangleEdgePair(crack_front_quadrangle_edges[hash_table_ind[j]],
                                                     crack_front_quadrangle_edges[hash_table_ind[k]],
                                                     -1.0 + EFFECTIVE_ZERO))
                {
                    crack_front_quadrangle_edges[hash_table_ind[j]]
                        = -crack_front_quadrangle_edges[hash_table_ind[j]] - 1;
                    crack_front_quadrangle_edges[hash_table_ind[k]]
                        = -crack_front_quadrangle_edges[hash_table_ind[k]] - 1;
                    break;
                }
            }
        }
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (crack_front_quadrangle_edges[hash_table_ind[j]] >= 0) continue;

            if (isZeroLengthQuadrangleEdge(-crack_front_quadrangle_edges[hash_table_ind[j]] - 1))
                crack_front_quadrangle_edges[hash_table_ind[j]]
                    = -crack_front_quadrangle_edges[hash_table_ind[j]] - 1;
        }
    offset = 0;
    for (i = 0; i < num_crack_front_quadrangle_edges; i++)
        if (crack_front_quadrangle_edges[i] < 0)
            crack_front_quadrangle_edges[offset++]
                = -crack_front_quadrangle_edges[i] - 1;
    num_crack_front_quadrangle_edges = offset;

    /* Destroy hahses and hash table */
    free(hash_table_ptr);
    free(hash_table_ind);
    free(hashes);
}
