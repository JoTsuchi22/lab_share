/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <autogl.h>

#include "edge.h"

#include "constant.h"
#include "scalar.h"

extern int division_number;

extern int **triangles;
extern int **quadrangles;
extern int num_visible_triangles;
extern int num_visible_quadrangles;
extern int *visible_triangles;
extern int *visible_quadrangles;
extern int num_triangle_edges;
extern int num_quadrangle_edges;
extern int *triangle_edges;
extern int *quadrangle_edges;
extern int num_nodes;
extern int num_dimensions;
extern double **nodes;

static void calculateTriangleNormalVector(double vector[3],
                                          int triangle_id);
static void calculateQuadrangleNormalVector(double vector[3],
                                            int quadrangle_id);

/*
 * Init or reinit edges
 */
void initOrReinitEdges(void)
{
    static int previous_num_triangle_edges = -1;
    static int previous_num_quadrangle_edges = -1;

    int i;

    /* Get number of edges */
    num_triangle_edges   = 3 * num_visible_triangles;
    num_quadrangle_edges = 4 * num_visible_quadrangles;

    /* Allocate edges */
    if (num_triangle_edges > previous_num_triangle_edges)
    {
        triangle_edges
            = (int *)realloc(triangle_edges,
                             sizeof(int) * num_triangle_edges);

        previous_num_triangle_edges = num_triangle_edges;
    }
    if (num_quadrangle_edges > previous_num_quadrangle_edges)
    {
        quadrangle_edges
            = (int *)realloc(quadrangle_edges,
                             sizeof(int) * num_quadrangle_edges);

        previous_num_quadrangle_edges = num_quadrangle_edges;
    }

    /* Generate edges */
    for (i = 0; i < num_visible_triangles; i++)
    {
        triangle_edges[3 * i    ] = 3 * visible_triangles[i];
        triangle_edges[3 * i + 1] = 3 * visible_triangles[i] + 1;
        triangle_edges[3 * i + 2] = 3 * visible_triangles[i] + 2;
    }
    for (i = 0; i < num_visible_quadrangles; i++)
    {
        quadrangle_edges[4 * i    ] = 4 * visible_quadrangles[i];
        quadrangle_edges[4 * i + 1] = 4 * visible_quadrangles[i] + 1;
        quadrangle_edges[4 * i + 2] = 4 * visible_quadrangles[i] + 2;
        quadrangle_edges[4 * i + 3] = 4 * visible_quadrangles[i] + 3;
    }

    /* Remove inner edges */
    removeInnerTriangleEdges();
    removeInnerQuadrangleEdges();
}

/*
 * Destroy edges
 */
void destroyEdges(void)
{
    free(quadrangle_edges);
    free(triangle_edges);
}

/*
 * Remove inner triangle edges
 */
void removeInnerTriangleEdges(void)
{
    int *hashes;
    int max_hash;
    double hash_factor;

    int *hash_table_ind;
    int *hash_table_ptr;
    int hash_table_size = num_triangle_edges;

    int offset;
    int i, j, k;

    if (num_triangle_edges == 0) return;

    /* Allocate hashes and hash table */
    hashes = (int *)malloc(sizeof(int) * num_triangle_edges);
    hash_table_ind = (int *)malloc(sizeof(int) * num_triangle_edges);
    hash_table_ptr = (int *)malloc(sizeof(int) * (hash_table_size + 1));

    /* Calculate hashes */
    max_hash = 1;
    for (i = 0; i < num_triangle_edges; i++)
    {
        int triangle_id      = triangle_edges[i] / 3;
        int triangle_edge_id = triangle_edges[i] % 3;

        int node0_id = triangles[triangle_id][triangle_edge_id];
        int node1_id = triangles[triangle_id][(triangle_edge_id + 1) % 3];

        hashes[i]
            = node0_id
            + node1_id;

        if (max_hash < hashes[i])
            max_hash = hashes[i];
    }
    hash_factor = (double)(hash_table_size - 1) / (double)max_hash;
    for (i = 0; i < num_triangle_edges; i++)
        hashes[i] = (int)(hash_factor * (double)hashes[i]);

    /* Generate hash table */
    for (i = 0; i < hash_table_size + 1; i++)
        hash_table_ptr[i] = 0;
    for (i = 0; i < num_triangle_edges; i++)
        hash_table_ptr[hashes[i] + 1]++;
    for (i = 0; i < hash_table_size; i++)
        hash_table_ptr[i + 1] += hash_table_ptr[i];
    for (i = 0; i < num_triangle_edges; i++)
        hash_table_ind[i] = -1;
    for (i = 0; i < num_triangle_edges; i++)
        for (j = hash_table_ptr[hashes[i]];
             j < hash_table_ptr[hashes[i] + 1];
             j++)
            if (hash_table_ind[j] == -1)
            {
                hash_table_ind[j] = i;
                break;
            }

    /* Find visible triangle edges */
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (triangle_edges[hash_table_ind[j]] < 0) continue;

            for (k = j + 1; k < hash_table_ptr[i + 1]; k++)
            {
                if (triangle_edges[hash_table_ind[k]] < 0) continue;

                if (isSameTriangleEdgePair(triangle_edges[hash_table_ind[j]],
                                           triangle_edges[hash_table_ind[k]]))
                {
                    triangle_edges[hash_table_ind[j]]
                        = -triangle_edges[hash_table_ind[j]] - 1;
                    triangle_edges[hash_table_ind[k]]
                        = -triangle_edges[hash_table_ind[k]] - 1;
                    break;
                }
            }
        }
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (triangle_edges[hash_table_ind[j]] >= 0) continue;

            for (k = j + 1; k < hash_table_ptr[i + 1]; k++)
            {
                if (triangle_edges[hash_table_ind[k]] >= 0) continue;

                if (isSameTriangleEdgePair(-triangle_edges[hash_table_ind[j]] - 1,
                                           -triangle_edges[hash_table_ind[k]] - 1)
                    && isNonSmoothTriangleEdgePair(-triangle_edges[hash_table_ind[j]] - 1,
                                                   -triangle_edges[hash_table_ind[k]] - 1,
                                                   cos(0.5 * M_PI / (double)division_number)))
                {
                    triangle_edges[hash_table_ind[j]]
                        = -triangle_edges[hash_table_ind[j]] - 1;
                    triangle_edges[hash_table_ind[k]]
                        = -triangle_edges[hash_table_ind[k]] - 1;
                    break;
                }
            }
        }
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (triangle_edges[hash_table_ind[j]] < 0) continue;

            if (isZeroLengthTriangleEdge(triangle_edges[hash_table_ind[j]]))
                triangle_edges[hash_table_ind[j]]
                    = -triangle_edges[hash_table_ind[j]] - 1;
        }
    offset = 0;
    for (i = 0; i < num_triangle_edges; i++)
        if (triangle_edges[i] >= 0)
            triangle_edges[offset++] = triangle_edges[i];
    num_triangle_edges = offset;

    /* Destroy hahses and hash table */
    free(hash_table_ptr);
    free(hash_table_ind);
    free(hashes);
}

/*
 * Remove inner quadrangle edges
 */
void removeInnerQuadrangleEdges(void)
{
    int *hashes;
    int max_hash;
    double hash_factor;

    int *hash_table_ind;
    int *hash_table_ptr;
    int hash_table_size = num_quadrangle_edges;

    int offset;
    int i, j, k;

    if (num_quadrangle_edges == 0) return;

    /* Allocate hashes and hash table */
    hashes = (int *)malloc(sizeof(int) * num_quadrangle_edges);
    hash_table_ind = (int *)malloc(sizeof(int) * num_quadrangle_edges);
    hash_table_ptr = (int *)malloc(sizeof(int) * (hash_table_size + 1));

    /* Calculate hashes */
    max_hash = 1;
    for (i = 0; i < num_quadrangle_edges; i++)
    {
        int quadrangle_id      = quadrangle_edges[i] / 4;
        int quadrangle_edge_id = quadrangle_edges[i] % 4;

        int node0_id = quadrangles[quadrangle_id][quadrangle_edge_id];
        int node1_id = quadrangles[quadrangle_id][(quadrangle_edge_id + 1) % 4];

        hashes[i]
            = node0_id
            + node1_id;

        if (max_hash < hashes[i])
            max_hash = hashes[i];
    }
    hash_factor = (double)(hash_table_size - 1) / (double)max_hash;
    for (i = 0; i < num_quadrangle_edges; i++)
        hashes[i] = (int)(hash_factor * (double)hashes[i]);

    /* Generate hash table */
    for (i = 0; i < hash_table_size + 1; i++)
        hash_table_ptr[i] = 0;
    for (i = 0; i < num_quadrangle_edges; i++)
        hash_table_ptr[hashes[i] + 1]++;
    for (i = 0; i < hash_table_size; i++)
        hash_table_ptr[i + 1] += hash_table_ptr[i];
    for (i = 0; i < num_quadrangle_edges; i++)
        hash_table_ind[i] = -1;
    for (i = 0; i < num_quadrangle_edges; i++)
        for (j = hash_table_ptr[hashes[i]];
             j < hash_table_ptr[hashes[i] + 1];
             j++)
            if (hash_table_ind[j] == -1)
            {
                hash_table_ind[j] = i;
                break;
            }

    /* Find visible quadrangle edges */
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (quadrangle_edges[hash_table_ind[j]] < 0) continue;

            for (k = j + 1; k < hash_table_ptr[i + 1]; k++)
            {
                if (quadrangle_edges[hash_table_ind[k]] < 0) continue;

                if (isSameQuadrangleEdgePair(quadrangle_edges[hash_table_ind[j]],
                                             quadrangle_edges[hash_table_ind[k]]))
                {
                    quadrangle_edges[hash_table_ind[j]]
                        = -quadrangle_edges[hash_table_ind[j]] - 1;
                    quadrangle_edges[hash_table_ind[k]]
                        = -quadrangle_edges[hash_table_ind[k]] - 1;
                    break;
                }
            }
        }
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (quadrangle_edges[hash_table_ind[j]] >= 0) continue;

            for (k = j + 1; k < hash_table_ptr[i + 1]; k++)
            {
                if (quadrangle_edges[hash_table_ind[k]] >= 0) continue;

                if (isSameQuadrangleEdgePair(-quadrangle_edges[hash_table_ind[j]] - 1,
                                             -quadrangle_edges[hash_table_ind[k]] - 1)
                    && isNonSmoothQuadrangleEdgePair(-quadrangle_edges[hash_table_ind[j]] - 1,
                                                     -quadrangle_edges[hash_table_ind[k]] - 1,
                                                     cos(0.5 * M_PI / (double)division_number)))
                {
                    quadrangle_edges[hash_table_ind[j]]
                        = -quadrangle_edges[hash_table_ind[j]] - 1;
                    quadrangle_edges[hash_table_ind[k]]
                        = -quadrangle_edges[hash_table_ind[k]] - 1;
                    break;
                }
            }
        }
    for (i = 0; i < hash_table_size; i++)
        for (j = hash_table_ptr[i]; j < hash_table_ptr[i + 1]; j++)
        {
            if (quadrangle_edges[hash_table_ind[j]] < 0) continue;

            if (isZeroLengthQuadrangleEdge(quadrangle_edges[hash_table_ind[j]]))
                quadrangle_edges[hash_table_ind[j]]
                    = -quadrangle_edges[hash_table_ind[j]] - 1;
        }
    offset = 0;
    for (i = 0; i < num_quadrangle_edges; i++)
        if (quadrangle_edges[i] >= 0)
            quadrangle_edges[offset++] = quadrangle_edges[i];
    num_quadrangle_edges = offset;

    /* Destroy hahses and hash table */
    free(hash_table_ptr);
    free(hash_table_ind);
    free(hashes);
}

/*
 * Is same triangle edge pair?
 */
int isSameTriangleEdgePair(int edge_a_id, int edge_b_id)
{
    int triangle_a_id = edge_a_id / 3;
    int triangle_b_id = edge_b_id / 3;

    int triangle_edge_a_id = edge_a_id % 3;
    int triangle_edge_b_id = edge_b_id % 3;

    int triangle_node_a0_id = triangles[triangle_a_id][triangle_edge_a_id];
    int triangle_node_a1_id = triangles[triangle_a_id][(triangle_edge_a_id + 1) % 3];
    int triangle_node_b0_id = triangles[triangle_b_id][triangle_edge_b_id];
    int triangle_node_b1_id = triangles[triangle_b_id][(triangle_edge_b_id + 1) % 3];

    sortTwoIntegers(&triangle_node_a0_id, &triangle_node_a1_id);
    sortTwoIntegers(&triangle_node_b0_id, &triangle_node_b1_id);

    return
        triangle_node_a0_id == triangle_node_b0_id
        && triangle_node_a1_id == triangle_node_b1_id;
}

/*
 * Is same quadrangle edge pair?
 */
int isSameQuadrangleEdgePair(int edge_a_id, int edge_b_id)
{
    int quadrangle_a_id = edge_a_id / 4;
    int quadrangle_b_id = edge_b_id / 4;

    int quadrangle_edge_a_id = edge_a_id % 4;
    int quadrangle_edge_b_id = edge_b_id % 4;

    int quadrangle_node_a0_id = quadrangles[quadrangle_a_id][quadrangle_edge_a_id];
    int quadrangle_node_a1_id = quadrangles[quadrangle_a_id][(quadrangle_edge_a_id + 1) % 4];
    int quadrangle_node_b0_id = quadrangles[quadrangle_b_id][quadrangle_edge_b_id];
    int quadrangle_node_b1_id = quadrangles[quadrangle_b_id][(quadrangle_edge_b_id + 1) % 4];

    sortTwoIntegers(&quadrangle_node_a0_id, &quadrangle_node_a1_id);
    sortTwoIntegers(&quadrangle_node_b0_id, &quadrangle_node_b1_id);

    return
        quadrangle_node_a0_id == quadrangle_node_b0_id
        && quadrangle_node_a1_id == quadrangle_node_b1_id;
}

/*
 * Is non-smooth triangle edge pair?
 */
int isNonSmoothTriangleEdgePair(int edge_a_id, int edge_b_id,
                                double threshold)
{
    double vector_a[3], vector_b[3];

    int triangle_a_id = edge_a_id / 3;
    int triangle_b_id = edge_b_id / 3;

    calculateTriangleNormalVector(vector_a, triangle_a_id);
    calculateTriangleNormalVector(vector_b, triangle_b_id);

    return
        AutoGL_GetScalarProduct3D(vector_a[0], vector_a[1], vector_a[2],
                                  vector_b[0], vector_b[1], vector_b[2])
        <= threshold;
}

/*
 * Is non-smooth quadrangle edge pair?
 */
int isNonSmoothQuadrangleEdgePair(int edge_a_id, int edge_b_id,
                                  double threshold)
{
    double vector_a[3], vector_b[3];

    int quadrangle_a_id = edge_a_id / 4;
    int quadrangle_b_id = edge_b_id / 4;

    calculateQuadrangleNormalVector(vector_a, quadrangle_a_id);
    calculateQuadrangleNormalVector(vector_b, quadrangle_b_id);

    return
        AutoGL_GetScalarProduct3D(vector_a[0], vector_a[1], vector_a[2],
                                  vector_b[0], vector_b[1], vector_b[2])
        <= threshold;
}

/*
 * Calculate triangle normal vector
 */
static void calculateTriangleNormalVector(double vector[3],
                                          int triangle_id)
{
    double a_x = nodes[triangles[triangle_id][0]][0];
    double a_y = nodes[triangles[triangle_id][0]][1];
    double a_z = nodes[triangles[triangle_id][0]][2];

    double b_x = nodes[triangles[triangle_id][1]][0];
    double b_y = nodes[triangles[triangle_id][1]][1];
    double b_z = nodes[triangles[triangle_id][1]][2];

    double c_x = nodes[triangles[triangle_id][2]][0];
    double c_y = nodes[triangles[triangle_id][2]][1];
    double c_z = nodes[triangles[triangle_id][2]][2];

    AutoGL_GetVectorProduct3D(&vector[0], &vector[1], &vector[2],
                              b_x - a_x, b_y - a_y, b_z - a_z,
                              c_x - a_x, c_y - a_y, c_z - a_z);
    AutoGL_NormalizeVector3D(&vector[0], &vector[1], &vector[2],
                             vector[0], vector[1], vector[2]);
}

/*
 * Calculate quadrangle normal vector
 */
static void calculateQuadrangleNormalVector(double vector[3],
                                            int quadrangle_id)
{
    double a_x = nodes[quadrangles[quadrangle_id][0]][0];
    double a_y = nodes[quadrangles[quadrangle_id][0]][1];
    double a_z = nodes[quadrangles[quadrangle_id][0]][2];

    double b_x = nodes[quadrangles[quadrangle_id][1]][0];
    double b_y = nodes[quadrangles[quadrangle_id][1]][1];
    double b_z = nodes[quadrangles[quadrangle_id][1]][2];

    double c_x = nodes[quadrangles[quadrangle_id][2]][0];
    double c_y = nodes[quadrangles[quadrangle_id][2]][1];
    double c_z = nodes[quadrangles[quadrangle_id][2]][2];

    double d_x = nodes[quadrangles[quadrangle_id][3]][0];
    double d_y = nodes[quadrangles[quadrangle_id][3]][1];
    double d_z = nodes[quadrangles[quadrangle_id][3]][2];

    double temp_x, temp_y, temp_z;

    AutoGL_GetVectorProduct3D(&vector[0], &vector[1], &vector[2],
                              b_x - a_x, b_y - a_y, b_z - a_z,
                              c_x - a_x, c_y - a_y, c_z - a_z);
    AutoGL_GetVectorProduct3D(&temp_x, &temp_y, &temp_z,
                              c_x - a_x, c_y - a_y, c_z - a_z,
                              d_x - a_x, d_y - a_y, d_z - a_z);
    AutoGL_AddVector3D(&vector[0], &vector[1], &vector[2],
                       vector[0], vector[1], vector[2],
                       temp_x, temp_y, temp_z);
    AutoGL_NormalizeVector3D(&vector[0], &vector[1], &vector[2],
                             vector[0], vector[1], vector[2]);
}

/*
 * Is zero-length triangle edge?
 */
int isZeroLengthTriangleEdge(int edge_id)
{
    int triangle_id      = edge_id / 3;
    int triangle_edge_id = edge_id % 3;

    int triangle_node0_id = triangles[triangle_id][triangle_edge_id];
    int triangle_node1_id = triangles[triangle_id][(triangle_edge_id + 1) % 3];

    return triangle_node0_id == triangle_node1_id;
}

/*
 * Is zero-length quadrangle edge?
 */
int isZeroLengthQuadrangleEdge(int edge_id)
{
    int quadrangle_id      = edge_id / 4;
    int quadrangle_edge_id = edge_id % 4;

    int quadrangle_node0_id = quadrangles[quadrangle_id][quadrangle_edge_id];
    int quadrangle_node1_id = quadrangles[quadrangle_id][(quadrangle_edge_id + 1) % 4];

    return quadrangle_node0_id == quadrangle_node1_id;
}
