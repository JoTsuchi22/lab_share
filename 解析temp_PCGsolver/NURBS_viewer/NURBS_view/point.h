/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#ifndef POINT_H
#define POINT_H

/*
 * Init or reinit points
 */
void initOrReinitPoints(void);

/*
 * Init or reinit points fields
 */
void initOrReinitPointsFields(void);

/*
 * Destroy points
 */
void destroyPoints(void);

/*
 * Destroy points fields
 */
void destroyPointsFields(void);

/*
 * Set points field
 */
void setPointsField(void);

/*
 * Set points field range
 */
void setPointsFieldRange(void);

#endif  /* POINT_H */
