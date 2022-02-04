/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#ifndef MESH_H
#define MESH_H

/*
 * Init or reinit mesh
 */
void initOrReinitMesh(void);

/*
 * Destroy mesh
 */
void destroyMesh(void);

/*
 * Reinit mesh from visible and invisible elements
 */
void reinitMeshFromVisibleAndInvisibleElements(void);

#endif  /* MESH_H */
