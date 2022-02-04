/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#ifndef FACE_H
#define FACE_H

/*
 * Init or reinit faces
 */
void initOrReinitFaces(void);

/*
 * Destroy faces
 */
void destroyFaces(void);

/*
 * Reinit faces from visible elements
 */
void reinitFacesFromVisibleElements(void);

/*
 * Remove inner triangles
 */
void removeInnerTriangles(void);

/*
 * Remove inner quadrangles
 */
void removeInnerQuadrangles(void);

#endif  /* FACE_H */
