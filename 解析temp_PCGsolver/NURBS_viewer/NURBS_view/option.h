/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#ifndef OPTION_H
#define OPTION_H

/*
 * Set options
 */
void setOptions(int argc, char *argv[]);

/*
 * Set options from arguments
 */
void setOptionsFromArguments(int argc, char *argv[]);

/*
 * Display help
 */
void displayHelp(void);

/*
 * Input restart file
 */
void inputRestartFile(void);

/*
 * Output restart file
 */
void outputRestartFile(void);

#endif  /* OPTION_H */
