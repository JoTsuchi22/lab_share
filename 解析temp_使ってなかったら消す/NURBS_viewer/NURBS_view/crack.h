/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#ifndef CRACK_H
#define CRACK_H

/*
 * Init or reinit crack faces
 */
void initOrReinitCrackFaces(void);

/*
 * Destroy crack faces
 */
void destroyCrackFaces(void);

/*
 * Remove non-crack triangles
 */
void removeNonCrackTriangles(void);

/*
 * Remove non-crack quadrangles
 */
void removeNonCrackQuadrangles(void);

/*
 * Init or reinit crack front edges
 */
void initOrReinitCrackFrontEdges(void);

/*
 * Destroy crack front edges
 */
void destroyCrackFrontEdges(void);

/*
 * Remove non-crack-front triangle edges
 */
void removeNonCrackFrontTriangleEdges(void);

/*
 * Remove non-crack-front quadrangle edges
 */
void removeNonCrackFrontQuadrangleEdges(void);

#endif  /* CRACK_H */
