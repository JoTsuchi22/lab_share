/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#ifndef EDGE_H
#define EDGE_H

/*
 * Init or reinit edges
 */
void initOrReinitEdges(void);

/*
 * Destroy edges
 */
void destroyEdges(void);

/*
 * Remove inner triangle edges
 */
void removeInnerTriangleEdges(void);

/*
 * Remove inner quadrangle edges
 */
void removeInnerQuadrangleEdges(void);

/*
 * Is same triangle edge pair?
 */
int isSameTriangleEdgePair(int edge_a_id, int edge_b_id);

/*
 * Is same quadrangle edge pair?
 */
int isSameQuadrangleEdgePair(int edge_a_id, int edge_b_id);

/*
 * Is non-smooth triangle edge pair?
 */
int isNonSmoothTriangleEdgePair(int edge_a_id, int edge_b_id,
                                double threshold);

/*
 * Is non-smooth quadrangle edge pair?
 */
int isNonSmoothQuadrangleEdgePair(int edge_a_id, int edge_b_id,
                                  double threshold);

/*
 * Is zero-length triangle edge?
 */
int isZeroLengthTriangleEdge(int edge_id);

/*
 * Is zero-length quadrangle edge?
 */
int isZeroLengthQuadrangleEdge(int edge_id);

#endif  /* EDGE_H */
