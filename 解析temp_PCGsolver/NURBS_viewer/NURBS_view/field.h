/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#ifndef FIELD_H
#define FIELD_H

/*
 * Init fields
 */
void initFields(void);

/*
 * Reinit fields
 */
void reinitFields(void);

/*
 * Destroy fields
 */
void destroyFields(void);

/*
 * Set contour field
 */
void setContourField(void);

/*
 * Set contour field range
 */
void setContourFieldRange(void);

/*
 * Set displacement field
 */
void setDisplacementField(void);

#endif  /* FIELD_H */
