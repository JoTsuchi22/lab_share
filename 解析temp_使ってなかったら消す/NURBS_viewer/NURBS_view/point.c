/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "point.h"

#include "constant.h"

extern char filename_points[];
extern char filename_points_fields[];

extern int current_time_step;

extern int num_dimensions;

extern int num_points;
extern double **points;
extern int num_points_fields;
extern double **points_fields;
extern char **points_field_labels;

extern int points_type;
extern char points_contour_bar_label[];
extern char points_contour_bar_label2[];
extern int auto_points_field_range_flag;
extern double points_field_max;
extern double points_field_min;
extern int points_field_id;

/*
 * Init or reinit points
 */
void initOrReinitPoints(void)
{
    FILE *fp;
    fpos_t position;

    char filename[256];
    char buffer[FILE_INPUT_BUFFER_SIZE];

    static char previous_filename[256] = "";

    void *flag;

    int num_point_dimensions;
    double point_coords[3];

    static int previous_num_points = -1;

    int i;

    if (strlen(filename_points) == 0) return;

    /* Open file */
    sprintf(filename, filename_points,
            current_time_step, current_time_step, current_time_step);
    if (strcmp(filename, previous_filename) == 0) return;
    strcpy(previous_filename, filename);
    printf("Input %s\n", filename);
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    /* Read number of points */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_points) != 1)
    {
        fprintf(stderr, "Error: failed to read number of points in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    printf("Number of points: %d\n", num_points);

    /* Read 1st point */
    fgetpos(fp, &position);
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
    {
        fprintf(stderr, "Error: failed to read point 0 in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }
    fsetpos(fp, &position);
    num_point_dimensions = sscanf(buffer,
                                  "%lf%lf%lf",
                                  &point_coords[0],
                                  &point_coords[1],
                                  &point_coords[2]);
    if (num_point_dimensions != num_dimensions)
    {
        fprintf(stderr, "Error: unsupported number of point dimensions (%d)\n",
                num_point_dimensions);
        exit(EXIT_FAILURE);
    }

    /* Allocate points */
    if (num_points > previous_num_points)
    {
        flag = (void *)points;
        points
            = (double **)realloc(points,
                                 sizeof(double *) * num_points);
        if (num_points)
        {
            if (flag == NULL) points[0] = NULL;
            points[0]
                = (double *)realloc(points[0],
                                    sizeof(double) * 3 * num_points);
        }
        for (i = 0; i < num_points; i++)
            points[i]
                = points[0] + 3 * i;

        previous_num_points = num_points;
    }

    /* Read points */
    for (i = 0; i < num_points; i++)
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer,
                      "%lf%lf%lf",
                      &points[i][0],
                      &points[i][1],
                      &points[i][2]) != num_dimensions)
        {
            fprintf(stderr, "Error: failed to read point %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

    /* Close */
    fclose(fp);

    /* Zero-fill coordinate Z if 2D */
    if (num_dimensions == 2)
        for (i = 0; i < num_points; i++)
            points[i][2] = 0.0;
}

/*
 * Init or reinit points fields
 */
void initOrReinitPointsFields(void)
{
    FILE *fp;
    fpos_t position;

    char filename[256];
    char buffer[FILE_INPUT_BUFFER_SIZE];

    static char previous_filename[256] = "";

    char *key, *value;

    char label[256];
    double fields[24];

    void *flag;

    static int previous_num_points = -1;
    static int previous_num_points_fields = -1;

    int num_items;
    int temp;
    int i, j;

    if (strlen(filename_points_fields) == 0) return;

    /* Open file */
    sprintf(filename, filename_points_fields,
            current_time_step, current_time_step, current_time_step);
    if (strcmp(filename, previous_filename) == 0) return;
    strcpy(previous_filename, filename);
    printf("Input %s\n", filename);
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    /* Read properties */
    strcpy(label, "NoLabel");
    num_items = 0;
    while (1)
    {
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to read properties in %s\n",
                    filename);
            exit(EXIT_FAILURE);
        }

        key   = strtok(buffer, "=\n");
        value = strtok(NULL, "=\n");
        if (key == NULL || value == NULL)
            break;

        if (strcmp(key, "label") == 0)
            strcpy(label, value);
        else if (strcmp(key, "num_items") == 0)
            num_items = atoi(value);
    }

    /* Check number of items */
    if (num_items != num_points)
    {
        fprintf(stderr, "Error: number of items is not number of points in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Read 1st fields */
    fgetpos(fp, &position);
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
    {
        fprintf(stderr, "Error: failed to read %s 0 in %s\n",
                label, filename);
        exit(EXIT_FAILURE);
    }
    fsetpos(fp, &position);
    num_points_fields
        = sscanf(buffer,
                 "%d:"
                 "%lf%lf%lf%lf%lf%lf"
                 "%lf%lf%lf%lf%lf%lf"
                 "%lf%lf%lf%lf%lf%lf"
                 "%lf%lf%lf%lf%lf%lf",
                 &temp,
                 &fields[0],
                 &fields[1],
                 &fields[2],
                 &fields[3],
                 &fields[4],
                 &fields[5],
                 &fields[6],
                 &fields[7],
                 &fields[8],
                 &fields[9],
                 &fields[10],
                 &fields[11],
                 &fields[12],
                 &fields[13],
                 &fields[14],
                 &fields[15],
                 &fields[16],
                 &fields[17],
                 &fields[18],
                 &fields[19],
                 &fields[20],
                 &fields[21],
                 &fields[22],
                 &fields[23]) - 1;

    /* Allocate points fields and points field labels */
    if (num_points > previous_num_points
        || num_points_fields > previous_num_points_fields)
    {
        flag = (void *)points_fields;
        points_fields
            = (double **)realloc(points_fields,
                                 sizeof(double *) * num_points);
        if (num_points_fields)
        {
            if (flag == NULL) points_fields[0] = NULL;
            points_fields[0]
                = (double *)realloc(points_fields[0],
                                    sizeof(double) * num_points_fields * num_points);
        }
        for (i = 0; i < num_points; i++)
            points_fields[i]
                = points_fields[0] + num_points_fields * i;
        flag = (void *)points_field_labels;
        points_field_labels
            = (char **)realloc(points_field_labels,
                               sizeof(char *) * num_points_fields);
        if (num_points_fields)
        {
            if (flag == NULL) points_field_labels[0] = NULL;
            points_field_labels[0]
                = (char *)realloc(points_field_labels[0],
                                  sizeof(char) * 256 * num_points_fields);
        }
        for (i = 0; i < num_points_fields; i++)
            points_field_labels[i]
                = points_field_labels[0] + 256 * i;

        previous_num_points = num_points;
        previous_num_points_fields = num_points_fields;
    }

    /* Get labels */
    if (num_points_fields == 1)
    {
        strcpy(points_field_labels[0],
               label);
    }
    else
    {
        for (i = 0; i < num_points_fields; i++)
            sprintf(points_field_labels[i],
                    "%s #%d", label, i + 1);
    }

    /* Read points fields */
    for (i = 0; i < num_points; i++)
    {
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer,
                      "%d:"
                      "%lf%lf%lf%lf%lf%lf"
                      "%lf%lf%lf%lf%lf%lf"
                      "%lf%lf%lf%lf%lf%lf"
                      "%lf%lf%lf%lf%lf%lf",
                      &temp,
                      &fields[0],
                      &fields[1],
                      &fields[2],
                      &fields[3],
                      &fields[4],
                      &fields[5],
                      &fields[6],
                      &fields[7],
                      &fields[8],
                      &fields[9],
                      &fields[10],
                      &fields[11],
                      &fields[12],
                      &fields[13],
                      &fields[14],
                      &fields[15],
                      &fields[16],
                      &fields[17],
                      &fields[18],
                      &fields[19],
                      &fields[20],
                      &fields[21],
                      &fields[22],
                      &fields[23])
            != num_points_fields + 1)
        {
            fprintf(stderr, "Error: failed to read %s %d in %s\n",
                    label, i, filename);
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < num_points_fields; j++)
            points_fields[i][j] = fields[j];
    }

    /* Close */
    fclose(fp);

    /* Set points field */
    if (points_field_id == -1
        && num_points_fields)
        points_field_id = 0;
    setPointsField();
    if (strlen(points_contour_bar_label2))
        strcpy(points_contour_bar_label,
               points_contour_bar_label2);

    /* Set points field range */
    if (auto_points_field_range_flag)
        setPointsFieldRange();
}

/*
 * Destroy points
 */
void destroyPoints(void)
{
    if (strlen(filename_points) == 0) return;

    if (num_points)
        free(points[0]);
    free(points);
}

/*
 * Destroy points fields
 */
void destroyPointsFields(void)
{
    if (strlen(filename_points_fields) == 0) return;

    if (num_points_fields)
        free(points_field_labels[0]);
    free(points_field_labels);

    if (num_points_fields)
        free(points_fields[0]);
    free(points_fields);
}

/*
 * Set points field
 */
void setPointsField(void)
{
    if (points_field_id >= 0
        && points_field_id < num_points_fields)
    {
        strcpy(points_contour_bar_label,
               points_field_labels[points_field_id]);
    }
    else
    {
        strcpy(points_contour_bar_label, "");
    }
}

/*
 * Set points field range
 */
void setPointsFieldRange(void)
{
    int i;

    if (points_field_id >= 0
        && points_field_id < num_points_fields)
    {
        points_field_max = -DBL_MAX;
        points_field_min = DBL_MAX;

        for (i = 0; i < num_points; i++)
        {
            if (points_field_max < points_fields[i][points_field_id])
                points_field_max = points_fields[i][points_field_id];
            if (points_field_min > points_fields[i][points_field_id])
                points_field_min = points_fields[i][points_field_id];
        }

        if (points_field_max <= points_field_min)
            points_field_max += 1.0;
    }
    else
    {
        points_field_max = 1.0;
        points_field_min = 0.0;
    }
}
