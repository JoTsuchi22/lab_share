/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#ifndef VIEW_H
#define VIEW_H

void calclateMinMax (void);

/*
 * Calculate view position from mesh
 */
void calculateViewPositionFromMesh(void);

/*
 * Align view directions
 */
void alignViewDirection(void);

/*
 * Align view up vector
 */
void alignViewUpVector(void);

/*
 * Step forward
 */
void stepForward(void);

/*
 * Step backward
 */
void stepBackward(void);

/*
 * Print view position
 */
void printViewPosition(void);

#endif  /* VIEW_H */
