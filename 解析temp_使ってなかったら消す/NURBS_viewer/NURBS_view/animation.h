/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#ifndef ANIMATION_H
#define ANIMATION_H

/*
 * Set last time steps
 */
void setTimeSteps(void);

/*
 * Play or pause
 */
void playOrPause(void);

/*
 * Increment time step
 */
void incrementTimeStep(void);

/*
 * Decrement time step
 */
void decrementTimeStep(void);

/*
 * Go to first time step
 */
void goToFirstTimeStep(void);

/*
 * Go to last time step
 */
void goToLastTimeStep(void);

/*
 * Save all view images
 */
void saveAllViewImages(void);

#endif  /* ANIMATION_H */
