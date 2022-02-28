/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <autogl.h>

#include "animation.h"

#include "mesh.h"
#include "field.h"
#include "point.h"
#include "image.h"

extern char *filename_mesh;
extern int num_filenames_fields;
extern char **filenames_fields;
extern char filename_visible_elements[];
extern char filename_invisible_elements[];
extern char filename_points[];
extern char filename_points_fields[];
extern char filename_view_image[];

extern int current_time_step;
extern int first_time_step;
extern int last_time_step;
extern int time_step_increment;
extern int time_step_loop_flag;
extern int slow_down_level;
extern char filename_all_view_images[];

static void incrementTimeStepAutomatically(void);
static void goToTimeStep(int time_step);

/*
 * Set time steps
 */
void setTimeSteps(void)
{
    FILE *fp;
    char filename[256];

    int last_time_step_mesh;
    int last_time_step_fields;
    int last_time_step_visible_elements;
    int last_time_step_invisible_elements;
    int last_time_step_points;
    int last_time_step_points_fields;

    int i, j;

    /* Check mesh files */
    last_time_step_mesh = INT_MAX;
    sprintf(filename, filename_mesh,
            first_time_step, first_time_step, first_time_step);
    if (strcmp(filename, filename_mesh))
    {
        for (i = first_time_step;; i += time_step_increment)
        {
            sprintf(filename, filename_mesh,
                    i, i, i);
            fp = fopen(filename, "r");
            if (fp == NULL) break;
            fclose(fp);
        }
        if (last_time_step_mesh > i - time_step_increment)
            last_time_step_mesh = i - time_step_increment;
    }

    /* Check fields files */
    last_time_step_fields = INT_MAX;
    for (i = 0; i < num_filenames_fields; i++)
    {
        sprintf(filename, filenames_fields[i],
                first_time_step, first_time_step, first_time_step);
        if (strcmp(filename, filenames_fields[i]))
        {
            for (j = first_time_step;; j += time_step_increment)
            {
                sprintf(filename, filenames_fields[i],
                        j, j, j);
                fp = fopen(filename, "r");
                if (fp == NULL) break;
                fclose(fp);
            }
            if (last_time_step_fields > j - time_step_increment)
                last_time_step_fields = j - time_step_increment;
        }
    }

    /* Check visible elements files */
    last_time_step_visible_elements = INT_MAX;
    if (strlen(filename_visible_elements))
    {
        sprintf(filename, filename_visible_elements,
                first_time_step, first_time_step, first_time_step);
        if (strcmp(filename, filename_visible_elements))
        {
            for (i = first_time_step;; i += time_step_increment)
            {
                sprintf(filename, filename_visible_elements,
                        i, i, i);
                fp = fopen(filename, "r");
                if (fp == NULL) break;
                fclose(fp);
            }
            if (last_time_step_visible_elements > i - time_step_increment)
                last_time_step_visible_elements = i - time_step_increment;
        }
    }

    /* Check invisible elements files */
    last_time_step_invisible_elements = INT_MAX;
    if (strlen(filename_invisible_elements))
    {
        sprintf(filename, filename_invisible_elements,
                first_time_step, first_time_step, first_time_step);
        if (strcmp(filename, filename_invisible_elements))
        {
            for (i = first_time_step;; i += time_step_increment)
            {
                sprintf(filename, filename_invisible_elements,
                        i, i, i);
                fp = fopen(filename, "r");
                if (fp == NULL) break;
                fclose(fp);
            }
            if (last_time_step_invisible_elements > i - time_step_increment)
                last_time_step_invisible_elements = i - time_step_increment;
        }
    }

    /* Check points files */
    last_time_step_points = INT_MAX;
    if (strlen(filename_points))
    {
        sprintf(filename, filename_points,
                first_time_step, first_time_step, first_time_step);
        if (strcmp(filename, filename_points))
        {
            for (i = first_time_step;; i += time_step_increment)
            {
                sprintf(filename, filename_points,
                        i, i, i);
                fp = fopen(filename, "r");
                if (fp == NULL) break;
                fclose(fp);
            }
            if (last_time_step_points > i - time_step_increment)
                last_time_step_points = i - time_step_increment;
        }
    }

    /* Check points fields files */
    last_time_step_points_fields = INT_MAX;
    if (strlen(filename_points_fields))
    {
        sprintf(filename, filename_points_fields,
                first_time_step, first_time_step, first_time_step);
        if (strcmp(filename, filename_points_fields))
        {
            for (i = first_time_step;; i += time_step_increment)
            {
                sprintf(filename, filename_points_fields,
                        i, i, i);
                fp = fopen(filename, "r");
                if (fp == NULL) break;
                fclose(fp);
            }
            if (last_time_step_points_fields > i - time_step_increment)
                last_time_step_points_fields = i - time_step_increment;
        }
    }

    /* Set last time step */
    last_time_step = INT_MAX;
    if (last_time_step > last_time_step_mesh)
        last_time_step = last_time_step_mesh;
    if (last_time_step > last_time_step_fields)
        last_time_step = last_time_step_fields;
    if (last_time_step > last_time_step_visible_elements)
        last_time_step = last_time_step_visible_elements;
    if (last_time_step > last_time_step_invisible_elements)
        last_time_step = last_time_step_invisible_elements;
    if (last_time_step > last_time_step_points)
        last_time_step = last_time_step_points;
    if (last_time_step > last_time_step_points_fields)
        last_time_step = last_time_step_points_fields;
    if (last_time_step == INT_MAX)
        last_time_step = first_time_step;

    /* Set current time step */
    current_time_step = first_time_step;

    if (first_time_step < last_time_step)
    {
        printf("First time step: %d\n", first_time_step);
        printf("Last time step: %d\n", last_time_step);
        printf("Time step increment: %d\n", time_step_increment);
        printf("Number of time steps: %d\n",
               (last_time_step - first_time_step) / time_step_increment + 1);
    }
}

/*
 * Play or pause
 */
void playOrPause(void)
{
    static int playing_flag = 0;

    if (playing_flag)
    {
        playing_flag = 0;
        AutoGL_SetIdleEventCallback(NULL);
    }
    else
    {
        playing_flag = 1;
        AutoGL_SetIdleEventCallback(incrementTimeStepAutomatically);
    }
}

/*
 * Increment time step automatically
 */
static void incrementTimeStepAutomatically(void)
{
    static int idle_count = 0;

    idle_count++;
    if (idle_count >= slow_down_level)
    {
        idle_count = 0;

        current_time_step += time_step_increment;
        if (current_time_step > last_time_step)
        {
            if (time_step_loop_flag)
            {
                current_time_step = first_time_step;
            }
            else
            {
                current_time_step = last_time_step;
                playOrPause();
            }
        }

        initOrReinitMesh();
        reinitMeshFromVisibleAndInvisibleElements();
        reinitFields();
        initOrReinitPoints();
        initOrReinitPointsFields();
    }

    AutoGL_DrawView();
}

/*
 * Increment time step
 */
void incrementTimeStep(void)
{
    goToTimeStep(current_time_step + time_step_increment);
}

/*
 * Decrement time step
 */
void decrementTimeStep(void)
{
    goToTimeStep(current_time_step - time_step_increment);
}

/*
 * Go to first time step
 */
void goToFirstTimeStep(void)
{
    goToTimeStep(first_time_step);
}

/*
 * Go to last time step
 */
void goToLastTimeStep(void)
{
    goToTimeStep(last_time_step);
}

/*
 * Go to time step
 */
static void goToTimeStep(int time_step)
{
    current_time_step = time_step;
    if (current_time_step < first_time_step)
    {
        current_time_step = first_time_step;
        return;
    }
    else if (current_time_step > last_time_step)
    {
        current_time_step = last_time_step;
        return;
    }

    initOrReinitMesh();
    reinitMeshFromVisibleAndInvisibleElements();
    reinitFields();
    initOrReinitPoints();
    initOrReinitPointsFields();

    AutoGL_DrawView();
}

/*
 * Save all view images
 */
void saveAllViewImages(void)
{
    char previous_filename_view_image[256];
    int previous_current_time_step;

    previous_current_time_step = current_time_step;
    strcpy(previous_filename_view_image,
           filename_view_image);

    for (current_time_step = first_time_step;
         current_time_step <= last_time_step;
         current_time_step += time_step_increment)
    {
        sprintf(filename_view_image, filename_all_view_images,
                current_time_step, current_time_step, current_time_step);

        initOrReinitMesh();
        reinitMeshFromVisibleAndInvisibleElements();
        reinitFields();
        initOrReinitPoints();
        initOrReinitPointsFields();

        AutoGL_DrawView();
        AutoGL_DrawView();    /* Somehow two times are necessary... */
        saveViewImage();
    }

    current_time_step = previous_current_time_step;
    strcpy(filename_view_image,
           previous_filename_view_image);

    initOrReinitMesh();
    reinitMeshFromVisibleAndInvisibleElements();
    reinitFields();
    initOrReinitPoints();
    initOrReinitPointsFields();

    AutoGL_DrawView();
}
