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

#ifdef WITH_ADVIO
#include <Adv/AdvDocument.h>
#endif

#include "field.h"

#include "constant.h"
#include "option.h"

extern int num_filenames_fields;
extern char **filenames_fields;

extern int current_time_step;
extern int first_time_step;
extern int last_time_step;
extern int time_step_increment;

extern int num_elements;
extern int num_segments;
extern int num_nodes;
extern int num_dimensions;

extern double *contour_field;
extern int num_nodal_fields;
extern double **nodal_fields;
extern char **nodal_field_labels;
extern int num_elemental_fields;
extern double **elemental_fields;
extern char **elemental_field_labels;
extern int num_integer_elemental_fields;
extern int **integer_elemental_fields;
extern double *casted_integer_elemental_field;
extern char **integer_elemental_field_labels;

extern double *displacement_fields[];
extern int num_displacement_fields;
extern int *displacement_field_offsets;
extern char **displacement_field_labels;

extern int integer_contour_field_mod;
extern char contour_bar_label[];
extern char contour_bar_label2[];
extern int auto_contour_field_range_flag;
extern double contour_field_max;
extern double contour_field_min;
extern int contour_field_id;
extern int displacement_field_id;

static double *zero_vector = NULL;

static void initOrReinitFields(void);
static void initOrReinitZeroVector(void);
static void initOrReinitFieldsFromAdvFile(void);
static void initOrReinitFieldsFromDatFiles(void);
static void initOrReinitFieldsFromMshFile(void);
static void initOrReinitFieldsFromEmshFile(void);
static void initOrReinitFieldsFromPcmFile(void);
static void initOrReinitFieldsFromVisibleElementsFile(void);

/*
 * Init fields
 */
void initFields(void)
{
    if (num_filenames_fields == 0) return;

    /* Init fields */
    initOrReinitFields();

    /* Set contour field */
    if (contour_field_id == -1)
    {
        if (num_nodal_fields)
            contour_field_id = 0;
        else if (num_elemental_fields)
            contour_field_id = num_nodal_fields;
        else if (num_integer_elemental_fields)
            contour_field_id = num_nodal_fields + num_elemental_fields;
        else
            contour_field_id = -1;
    }
    setContourField();
    if (strlen(contour_bar_label2))
        strcpy(contour_bar_label,
               contour_bar_label2);

    /* Set contour field range */
    if (auto_contour_field_range_flag)
        setContourFieldRange();

    /* Set displacement field */
    if (displacement_field_id == -1
        && num_displacement_fields)
        displacement_field_id = 0;
    setDisplacementField();
}

/*
 * Reinit fields
 */
void reinitFields(void)
{
    if (num_filenames_fields == 0) return;

    /* Reinit fields */
    initOrReinitFields();
}

/*
 * Init or reinit fields
 */
static void initOrReinitFields(void)
{
    char filename[256];

    static char previous_filename[256] = "";

    char *filename_suffix;
    int i;

    /* Check filename */
    sprintf(filename, filenames_fields[0],
            current_time_step, current_time_step, current_time_step);
    if (strcmp(filename, previous_filename) == 0) return;
    strcpy(previous_filename, filename);

    /* Get filename suffix */
    filename_suffix = strrchr(filenames_fields[0], '.');
    if (filename_suffix == NULL)
        filename_suffix = filenames_fields[0];

    /* Check number of filenames */
    if (strcmp(filename_suffix, ".dat") != 0
        && num_filenames_fields >= 2)
    {
        displayHelp();
        exit(EXIT_FAILURE);
    }
    if (strcmp(filename_suffix, ".dat") == 0)
        for (i = 1; i < num_filenames_fields; i++)
        {
            filename_suffix = strrchr(filenames_fields[i], '.');
            if (filename_suffix == NULL)
                filename_suffix = filenames_fields[i];

            if (strcmp(filename_suffix, ".dat"))
            {
                displayHelp();
                exit(EXIT_FAILURE);
            }
        }

    /* Init or reinit zero vector */
    initOrReinitZeroVector();

    /* Init or reinit fields */
    if (strcmp(filename_suffix, ".adv") == 0)
        initOrReinitFieldsFromAdvFile();
    else if (strcmp(filename_suffix, ".dat") == 0)
        initOrReinitFieldsFromDatFiles();
    else if (strcmp(filename_suffix, ".msh") == 0)
        initOrReinitFieldsFromMshFile();
    else if (strcmp(filename_suffix, ".emsh") == 0)
        initOrReinitFieldsFromEmshFile();
    else if (strcmp(filename_suffix, ".pcm") == 0)
        initOrReinitFieldsFromPcmFile();
    else
        initOrReinitFieldsFromVisibleElementsFile();
}

/*
 * Init or reinit zero vector
 */
static void initOrReinitZeroVector(void)
{
    static int previous_num = -1;

    int num;
    int i;

    num = (num_nodes > num_elements)? num_nodes: num_elements;
    if (num > previous_num)
    {
        zero_vector = (double *)realloc(zero_vector,
                                        sizeof(double) * num);
        for (i = 0; i < num; i++)
            zero_vector[i] = 0.0;

        previous_num = num;
    }
}

/*
 * Init or reinit fields from adv file
 */
static void initOrReinitFieldsFromAdvFile(void)
#ifdef WITH_ADVIO
{
    AdvDocFile *dfile;
    AdvDocument *doc;

    char filename[256];
    char label[256];
    char format[256];

    void *flag;

    static int previous_num_nodes = -1;
    static int previous_num_elements = -1;

    int previous_num_nodal_fields;
    int previous_num_displacement_fields;
    int previous_num_elemental_fields;
    int previous_num_integer_elemental_fields;
    int format_length;
    adv_off_t offset;
    int field_offset;
    int displacement_field_offset;
    int i, j, k;

    /* Open file */
    sprintf(filename, filenames_fields[0],
            current_time_step, current_time_step, current_time_step);
    printf("Input %s\n", filename);
    dfile = adv_dio_file_open(filename, "r");
    if (dfile == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Count nodal fields */
    previous_num_nodal_fields = num_nodal_fields;
    previous_num_displacement_fields = num_displacement_fields;
    num_nodal_fields = 0;
    num_displacement_fields = 0;
    for (i = 0;; i++)
    {
        int num_nodal_dofs;

        /* Open document */
        doc = adv_dio_open_nth(dfile, i);
        if (doc == NULL)
            break;

        /* Check properties */
        if (strcmp(adv_dio_get_property(doc, "content_type"), "FEGenericAttribute")
            || strcmp(adv_dio_get_property(doc, "fega_type"), "AllNodeVariable"))
            continue;

        /* Check format */
        strcpy(format, adv_dio_get_property(doc, "format"));
        format_length = (int)strlen(format);
        for (j = 0; j < format_length; j += 2)
            if (strncmp(&format[j], "f8", 2))
                break;
        if (j != format_length)
            continue;

        /* Get nodal DOFs */
        num_nodal_dofs
            = adv_format_get_size(adv_dio_get_property(doc, "format")) / 8;

        /* Count nodal fields */
        num_nodal_fields
            += num_nodal_dofs;

        /* Count displacement fields */
        if (num_nodal_dofs == num_dimensions)
            num_displacement_fields++;

        /* Close document */
        adv_dio_close(doc);
    }
    if (nodal_fields != NULL
        && num_nodal_fields != previous_num_nodal_fields)
    {
        fprintf(stderr, "Error: number of nodal fields changed\n");
        exit(EXIT_FAILURE);
    }
    if (displacement_fields[0] != NULL
        && num_displacement_fields != previous_num_displacement_fields)
    {
        fprintf(stderr, "Error: number of displacement fields changed\n");
        exit(EXIT_FAILURE);
    }

    /* Count elemental fields */
    previous_num_elemental_fields = num_elemental_fields;
    num_elemental_fields = 0;
    if (num_elements)
        for (i = 0;; i++)
        {
            int num_elemental_dofs;

            /* Open document */
            doc = adv_dio_open_nth(dfile, i);
            if (doc == NULL)
                break;

            /* Check properties */
            if (strcmp(adv_dio_get_property(doc, "content_type"), "FEGenericAttribute")
                || strcmp(adv_dio_get_property(doc, "fega_type"), "AllElementVariable"))
                continue;

            /* Check label */
            if (strstr(adv_dio_get_property(doc, "label"), "@IntegrationPoint"))
                continue;

            /* Check format */
            strcpy(format, adv_dio_get_property(doc, "format"));
            format_length = (int)strlen(format);
            for (j = 0; j < format_length; j += 2)
                if (strncmp(&format[j], "f8", 2))
                    break;
            if (j != format_length)
                continue;

            /* Get elemental DOFs */
            num_elemental_dofs
                = adv_format_get_size(adv_dio_get_property(doc, "format")) / 8;

            /* Count elemental fields */
            num_elemental_fields
                += num_elemental_dofs;

            /* Close document */
            adv_dio_close(doc);
        }
    if (elemental_fields != NULL
        && num_elemental_fields != previous_num_elemental_fields)
    {
        fprintf(stderr, "Error: number of elemental fields changed\n");
        exit(EXIT_FAILURE);
    }

    /* Count integer elemental fields */
    previous_num_integer_elemental_fields = num_integer_elemental_fields;
    num_integer_elemental_fields = 0;
    if (num_elements)
        for (i = 0;; i++)
        {
            /* Open document */
            doc = adv_dio_open_nth(dfile, i);
            if (doc == NULL)
                break;

            /* Check properties */
            if (strcmp(adv_dio_get_property(doc, "content_type"), "FEGenericAttribute")
                || strcmp(adv_dio_get_property(doc, "fega_type"), "AllElementVariable"))
                continue;

            /* Check label */
            if (strstr(adv_dio_get_property(doc, "label"), "@IntegrationPoint"))
                continue;

            /* Check format */
            strcpy(format, adv_dio_get_property(doc, "format"));
            format_length = (int)strlen(format);
            for (j = 0; j < format_length; j += 2)
                if (strncmp(&format[j], "i4", 2))
                    break;
            if (j != format_length)
                continue;

            /* Count elemental fields */
            num_integer_elemental_fields++;

            /* Close document */
            adv_dio_close(doc);
        }
    if (integer_elemental_fields != NULL
        && num_integer_elemental_fields != previous_num_integer_elemental_fields)
    {
        fprintf(stderr, "Error: number of integer elemental fields changed\n");
        exit(EXIT_FAILURE);
    }

    /* Allocate fields and field labels */
    if (num_nodes > previous_num_nodes)
    {
        flag = (void *)nodal_fields;
        nodal_fields
            = (double **)realloc(nodal_fields,
                                 sizeof(double *) * num_nodal_fields);
        if (num_nodal_fields)
        {
            if (flag == NULL) nodal_fields[0] = NULL;
            nodal_fields[0]
                = (double *)realloc(nodal_fields[0],
                                    sizeof(double) * num_nodes * num_nodal_fields);
        }
        for (i = 0; i < num_nodal_fields; i++)
            nodal_fields[i]
                = nodal_fields[0] + num_nodes * i;

        flag = (void *)nodal_field_labels;
        nodal_field_labels
            = (char **)realloc(nodal_field_labels,
                               sizeof(char *) * num_nodal_fields);
        if (num_nodal_fields)
        {
            if (flag == NULL) nodal_field_labels[0] = NULL;
            nodal_field_labels[0]
                = (char *)realloc(nodal_field_labels[0],
                                  sizeof(char) * 256 * num_nodal_fields);
        }
        for (i = 0; i < num_nodal_fields; i++)
            nodal_field_labels[i]
                = nodal_field_labels[0] + 256 * i;
    }
    if (num_elements > previous_num_elements)
    {
        flag = (void *)elemental_fields;
        elemental_fields
            = (double **)realloc(elemental_fields,
                                 sizeof(double *) * num_elemental_fields);
        if (num_elemental_fields)
        {
            if (flag == NULL) elemental_fields[0] = NULL;
            elemental_fields[0]
                = (double *)realloc(elemental_fields[0],
                                    sizeof(double) * num_elements * num_elemental_fields);
        }
        for (i = 0; i < num_elemental_fields; i++)
            elemental_fields[i]
                = elemental_fields[0] + num_elements * i;

        flag = (void *)elemental_field_labels;
        elemental_field_labels
            = (char **)realloc(elemental_field_labels,
                               sizeof(char *) * num_elemental_fields);
        if (num_elemental_fields)
        {
            if (flag == NULL) elemental_field_labels[0] = NULL;
            elemental_field_labels[0]
                = (char *)realloc(elemental_field_labels[0],
                                  sizeof(char) * 256 * num_elemental_fields);
        }
        for (i = 0; i < num_elemental_fields; i++)
            elemental_field_labels[i]
                = elemental_field_labels[0] + 256 * i;

        flag = (void *)integer_elemental_fields;
        integer_elemental_fields
            = (int **)realloc(integer_elemental_fields,
                              sizeof(int *) * num_integer_elemental_fields);
        if (num_integer_elemental_fields)
        {
            if (flag == NULL) integer_elemental_fields[0] = NULL;
            integer_elemental_fields[0]
                = (int *)realloc(integer_elemental_fields[0],
                                 sizeof(int) * num_elements * num_integer_elemental_fields);
        }
        for (i = 0; i < num_integer_elemental_fields; i++)
            integer_elemental_fields[i]
                = integer_elemental_fields[0] + num_elements * i;

        casted_integer_elemental_field
            = (double *)realloc(casted_integer_elemental_field,
                                sizeof(double) * num_elements);

        flag = (void *)integer_elemental_field_labels;
        integer_elemental_field_labels
            = (char **)realloc(integer_elemental_field_labels,
                               sizeof(char *) * num_integer_elemental_fields);
        if (num_integer_elemental_fields)
        {
            if (flag == NULL) integer_elemental_field_labels[0] = NULL;
            integer_elemental_field_labels[0]
                = (char *)realloc(integer_elemental_field_labels[0],
                                  sizeof(char) * 256 * num_integer_elemental_fields);
        }
        for (i = 0; i < num_integer_elemental_fields; i++)
            integer_elemental_field_labels[i]
                = integer_elemental_field_labels[0] + 256 * i;
    }

    /* Allocate displacement fields, offsets and labels */
    if (num_nodes > previous_num_nodes)
        if (num_displacement_fields)
        {
            displacement_field_offsets
                = (int *)realloc(displacement_field_offsets,
                                 sizeof(int) * num_displacement_fields);

            flag = (void *)displacement_field_labels;
            displacement_field_labels
                = (char **)realloc(displacement_field_labels,
                                   sizeof(char *) * num_displacement_fields);
            if (num_displacement_fields)
            {
                if (flag == NULL) displacement_field_labels[0] = NULL;
                displacement_field_labels[0]
                    = (char *)realloc(displacement_field_labels[0],
                                      sizeof(char) * 256 * num_displacement_fields);
            }
            for (i = 0; i < num_displacement_fields; i++)
                displacement_field_labels[i]
                    = displacement_field_labels[0] + 256 * i;
        }

    if (num_nodes > previous_num_nodes)
        previous_num_nodes = num_nodes;
    if (num_elements > previous_num_elements)
        previous_num_elements = num_elements;

    /* Read nodal fields */
    field_offset = 0;
    displacement_field_offset = 0;
    for (i = 0;; i++)
    {
        int num_nodal_dofs;

        /* Open document */
        doc = adv_dio_open_nth(dfile, i);
        if (doc == NULL)
            break;

        /* Check properties */
        if (strcmp(adv_dio_get_property(doc, "content_type"), "FEGenericAttribute")
            || strcmp(adv_dio_get_property(doc, "fega_type"), "AllNodeVariable"))
            continue;

        /* Check format */
        strcpy(format, adv_dio_get_property(doc, "format"));
        format_length = (int)strlen(format);
        for (j = 0; j < format_length; j += 2)
            if (strncmp(&format[j], "f8", 2))
                break;
        if (j != format_length)
            continue;

        /* Get nodal DOFs */
        num_nodal_dofs
            = adv_format_get_size(adv_dio_get_property(doc, "format")) / 8;

        /* Get labels */
        strcpy(label, adv_dio_get_property(doc, "label"));
        if (num_nodal_dofs == 1)
        {
            strcpy(nodal_field_labels[field_offset],
                   label);
        }
        else
        {
            for (j = 0; j < num_nodal_dofs; j++)
                sprintf(nodal_field_labels[field_offset + j],
                        "%s #%d", label, j + 1);
        }

        /* Get displacement field */
        if (num_nodal_dofs == num_dimensions)
        {
            displacement_field_offsets[displacement_field_offset]
                = field_offset;
            strcpy(displacement_field_labels[displacement_field_offset],
                   label);
            displacement_field_offset++;
        }

        /* Read data */
        offset = 0;
        for (j = 0; j < num_nodes; j++)
            for (k = 0; k < num_nodal_dofs; k++)
                offset
                    += adv_dio_read_float64(doc, offset,
                                            &nodal_fields[field_offset + k][j]);
        field_offset += num_nodal_dofs;

        /* Close document */
        adv_dio_close(doc);
    }

    /* Read elemental fields */
    field_offset = 0;
    if (num_elements)
        for (i = 0;; i++)
        {
            int num_elemental_dofs;

            /* Open document */
            doc = adv_dio_open_nth(dfile, i);
            if (doc == NULL)
                break;

            /* Check properties */
            if (strcmp(adv_dio_get_property(doc, "content_type"), "FEGenericAttribute")
                || strcmp(adv_dio_get_property(doc, "fega_type"), "AllElementVariable"))
                continue;

            /* Check label */
            if (strstr(adv_dio_get_property(doc, "label"), "@IntegrationPoint"))
                continue;

            /* Check format */
            strcpy(format, adv_dio_get_property(doc, "format"));
            format_length = (int)strlen(format);
            for (j = 0; j < format_length; j += 2)
                if (strncmp(&format[j], "f8", 2))
                    break;
            if (j != format_length)
                continue;

            /* Get elemental DOFs */
            num_elemental_dofs
                = adv_format_get_size(adv_dio_get_property(doc, "format")) / 8;

            /* Get labels */
            strcpy(label, adv_dio_get_property(doc, "label"));
            if (num_elemental_dofs == 1)
            {
                strcpy(elemental_field_labels[field_offset],
                       label);
            }
            else
            {
                for (j = 0; j < num_elemental_dofs; j++)
                    sprintf(elemental_field_labels[field_offset + j],
                            "%s #%d", label, j + 1);
            }

            /* Read data */
            offset = 0;
            for (j = 0; j < num_elements; j++)
                for (k = 0; k < num_elemental_dofs; k++)
                    offset
                        += adv_dio_read_float64(doc, offset,
                                                &elemental_fields[field_offset + k][j]);
            field_offset += num_elemental_dofs;

            /* Close document */
            adv_dio_close(doc);
        }

    /* Read integer elemental fields */
    field_offset = 0;
    if (num_elements)
        for (i = 0;; i++)
        {
            int num_integer_elemental_dofs;

            /* Open document */
            doc = adv_dio_open_nth(dfile, i);
            if (doc == NULL)
                break;

            /* Check properties */
            if (strcmp(adv_dio_get_property(doc, "content_type"), "FEGenericAttribute")
                || strcmp(adv_dio_get_property(doc, "fega_type"), "AllElementVariable"))
                continue;

            /* Check label */
            if (strstr(adv_dio_get_property(doc, "label"), "@IntegrationPoint"))
                continue;

            /* Check format */
            strcpy(format, adv_dio_get_property(doc, "format"));
            format_length = (int)strlen(format);
            for (j = 0; j < format_length; j += 2)
                if (strncmp(&format[j], "i4", 2))
                    break;
            if (j != format_length)
                continue;

            /* Get integer elemental DOFs */
            num_integer_elemental_dofs
                = adv_format_get_size(adv_dio_get_property(doc, "format")) / 4;

            /* Get labels */
            strcpy(label, adv_dio_get_property(doc, "label"));
            if (num_integer_elemental_dofs)
            {
                strcpy(integer_elemental_field_labels[field_offset],
                       label);
            }
            else
            {
                for (j = 0; j < num_integer_elemental_dofs; j++)
                    sprintf(integer_elemental_field_labels[field_offset + j],
                            "%s #%d", label, j + 1);
            }

            /* Read data */
            offset = 0;
            for (j = 0; j < num_elements; j++)
                offset
                    += adv_dio_read_int32(doc, offset,
                                          &integer_elemental_fields[field_offset][j]);
            field_offset++;

            /* Close document */
            adv_dio_close(doc);
        }

    /* Close file */
    adv_dio_file_close(dfile);
}
#else
{
    fprintf(stderr, "Error: adv file not supported in this build\n");
    exit(EXIT_FAILURE);
}
#endif

/*
 * Init or reinit fields from dat files
 */
static void initOrReinitFieldsFromDatFiles(void)
{
    FILE *fp;
    fpos_t position;

    char filename[256];
    char buffer[FILE_INPUT_BUFFER_SIZE];

    char *key, *value;

    char label[256];
    double fields[24];

    void *flag;

    static int previous_num_nodes = -1;
    static int previous_num_elements = -1;

    int previous_num_nodal_fields;
    int previous_num_elemental_fields;
    int previous_num_displacement_fields;
    int num_items;
    int num_nodal_or_elemental_dofs;
    int nodal_field_offset;
    int elemental_field_offset;
    int displacement_field_offset;
    int temp;
    int i, j, k;

    /* Count nodal and elemental fields */
    previous_num_nodal_fields = num_nodal_fields;
    previous_num_elemental_fields = num_elemental_fields;
    previous_num_displacement_fields = num_displacement_fields;
    num_nodal_fields = 0;
    num_elemental_fields = 0;
    num_displacement_fields = 0;
    for (i = 0; i < num_filenames_fields; i++)
    {
        /* Open file */
        sprintf(filename, filenames_fields[i],
                current_time_step, current_time_step, current_time_step);
        printf("Input %s\n", filename);
        fp = fopen(filename, "r");
        if (fp == NULL)
        {
            fprintf(stderr, "Error: failed to open file %s\n",
                    filename);
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
        if (num_items != num_nodes && num_items != num_elements)
        {
            fprintf(stderr, "Error: number of items is neither number of nodes nor number of elements in %s\n",
                    filename);
            exit(EXIT_FAILURE);
        }

        /* Read 1st fields */
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to read %s 0 in %s\n",
                    label, filename);
            exit(EXIT_FAILURE);
        }
        num_nodal_or_elemental_dofs
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

        /* Count nodal or elemental fields */
        if (num_items == num_nodes)
            num_nodal_fields += num_nodal_or_elemental_dofs;
        else if (num_items == num_elements)
            num_elemental_fields += num_nodal_or_elemental_dofs;

        if (num_nodes == num_elements)
            fprintf(stderr, "Warning: %s regarded as nodal field\n",
                    label);

        /* Count displacement fields */
        if (num_items == num_nodes
            && num_nodal_or_elemental_dofs == num_dimensions)
            num_displacement_fields++;

        /* Close */
        fclose(fp);
    }
    if (nodal_fields != NULL
        && num_nodal_fields != previous_num_nodal_fields)
    {
        fprintf(stderr, "Error: number of nodal fields changed\n");
        exit(EXIT_FAILURE);
    }
    if (elemental_fields != NULL
        && num_elemental_fields != previous_num_elemental_fields)
    {
        fprintf(stderr, "Error: number of elemental fields changed\n");
        exit(EXIT_FAILURE);
    }
    if (displacement_fields[0] != NULL
        && num_displacement_fields != previous_num_displacement_fields)
    {
        fprintf(stderr, "Error: number of displacement fields changed\n");
        exit(EXIT_FAILURE);
    }

    /* Allocate fields and field labels */
    if (num_nodes > previous_num_nodes)
    {
        flag = (void *)nodal_fields;
        nodal_fields
            = (double **)realloc(nodal_fields,
                                 sizeof(double *) * num_nodal_fields);
        if (num_nodal_fields)
        {
            if (flag == NULL) nodal_fields[0] = NULL;
            nodal_fields[0]
                = (double *)realloc(nodal_fields[0],
                                    sizeof(double) * num_nodes * num_nodal_fields);
        }
        for (i = 0; i < num_nodal_fields; i++)
            nodal_fields[i]
                = nodal_fields[0] + num_nodes * i;

        flag = (void *)nodal_field_labels;
        nodal_field_labels
            = (char **)realloc(nodal_field_labels,
                               sizeof(char *) * num_nodal_fields);
        if (num_nodal_fields)
        {
            if (flag == NULL) nodal_field_labels[0] = NULL;
            nodal_field_labels[0]
                = (char *)realloc(nodal_field_labels[0],
                                  sizeof(char) * 256 * num_nodal_fields);
        }
        for (i = 0; i < num_nodal_fields; i++)
            nodal_field_labels[i]
                = nodal_field_labels[0] + 256 * i;
    }
    if (num_elements > previous_num_elements)
    {
        flag = (void *)elemental_fields;
        elemental_fields
            = (double **)realloc(elemental_fields,
                                 sizeof(double *) * num_elemental_fields);
        if (num_elemental_fields)
        {
            if (flag == NULL) elemental_fields[0] = NULL;
            elemental_fields[0]
                = (double *)realloc(elemental_fields[0],
                                    sizeof(double) * num_elements * num_elemental_fields);
        }
        for (i = 0; i < num_elemental_fields; i++)
            elemental_fields[i]
                = elemental_fields[0] + num_elements * i;

        flag = (void *)elemental_field_labels;
        elemental_field_labels
            = (char **)realloc(elemental_field_labels,
                               sizeof(char *) * num_elemental_fields);
        if (num_elemental_fields)
        {
            if (flag == NULL) elemental_field_labels[0] = NULL;
            elemental_field_labels[0]
                = (char *)realloc(elemental_field_labels[0],
                                  sizeof(char) * 256 * num_elemental_fields);
        }
        for (i = 0; i < num_elemental_fields; i++)
            elemental_field_labels[i]
                = elemental_field_labels[0] + 256 * i;
    }

    /* Allocate displacement fields, offsets and labels */
    if (num_nodes > previous_num_nodes)
        if (num_displacement_fields)
        {
            displacement_field_offsets
                = (int *)realloc(displacement_field_offsets,
                                 sizeof(int) * num_displacement_fields);

            flag = (void *)displacement_field_labels;
            displacement_field_labels
                = (char **)realloc(displacement_field_labels,
                                   sizeof(char *) * num_displacement_fields);
            if (num_displacement_fields)
            {
                if (flag == NULL) displacement_field_labels[0] = NULL;
                displacement_field_labels[0]
                    = (char *)realloc(displacement_field_labels[0],
                                      sizeof(char) * 256 * num_displacement_fields);
            }
            for (i = 0; i < num_displacement_fields; i++)
                displacement_field_labels[i]
                    = displacement_field_labels[0] + 256 * i;
        }

    if (num_nodes > previous_num_nodes)
        previous_num_nodes = num_nodes;
    if (num_elements > previous_num_elements)
        previous_num_elements = num_elements;

    /* Read nodal and elemental fields */
    nodal_field_offset = 0;
    elemental_field_offset = 0;
    displacement_field_offset = 0;
    for (i = 0; i < num_filenames_fields; i++)
    {
        /* Open file */
        sprintf(filename, filenames_fields[i],
                current_time_step, current_time_step, current_time_step);
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
        if (num_items != num_nodes && num_items != num_elements)
        {
            fclose(fp);
            continue;
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
        num_nodal_or_elemental_dofs
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

        /* Get labels */
        if (num_items == num_nodes)
        {
            if (num_nodal_or_elemental_dofs == 1)
            {
                strcpy(nodal_field_labels[nodal_field_offset],
                       label);
            }
            else
            {
                for (j = 0; j < num_nodal_or_elemental_dofs; j++)
                    sprintf(nodal_field_labels[nodal_field_offset + j],
                            "%s #%d", label, j + 1);
            }
        }
        else if (num_items == num_elements)
        {
            if (num_nodal_or_elemental_dofs == 1)
            {
                strcpy(elemental_field_labels[elemental_field_offset],
                       label);
            }
            else
            {
                for (j = 0; j < num_nodal_or_elemental_dofs; j++)
                    sprintf(elemental_field_labels[elemental_field_offset + j],
                            "%s #%d", label, j + 1);
            }
        }

        /* Get displacement field */
        if (num_items == num_nodes
            && num_nodal_or_elemental_dofs == num_dimensions)
        {
            displacement_field_offsets[displacement_field_offset]
                = nodal_field_offset;
            strcpy(displacement_field_labels[displacement_field_offset],
                   label);
            displacement_field_offset++;
        }

        /* Read fields */
        if (num_items == num_nodes)
        {
            for (j = 0; j < num_nodes; j++)
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
                    != num_nodal_or_elemental_dofs + 1)
                {
                    fprintf(stderr, "Error: failed to read %s %d in %s\n",
                            label, j, filename);
                    exit(EXIT_FAILURE);
                }

                for (k = 0; k < num_nodal_or_elemental_dofs; k++)
                    nodal_fields[nodal_field_offset + k][j] = fields[k];
            }
        }
        else if (num_items == num_elements)
        {
            for (j = 0; j < num_elements; j++)
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
                    != num_nodal_or_elemental_dofs + 1)
                {
                    fprintf(stderr, "Error: failed to read %s %d in %s\n",
                            label, j, filename);
                    exit(EXIT_FAILURE);
                }

                for (k = 0; k < num_nodal_or_elemental_dofs; k++)
                    elemental_fields[elemental_field_offset + k][j]
                        = fields[k];
            }
        }

        /* Update offset */
        if (num_items == num_nodes)
            nodal_field_offset += num_nodal_or_elemental_dofs;
        else if (num_items == num_elements)
            elemental_field_offset += num_nodal_or_elemental_dofs;

        /* Close */
        fclose(fp);
    }
}

/*
 * Init or reinit fields from msh file
 */
static void initOrReinitFieldsFromMshFile(void)
{
    FILE *fp;

    char filename[256];
    char buffer[FILE_INPUT_BUFFER_SIZE];

    void *flag;

    static int previous_num_elements = -1;

    int num_items;
    int num_material_ids;
    int i, j;

    /* Open file */
    sprintf(filename, filenames_fields[0],
            current_time_step, current_time_step, current_time_step);
    printf("Input %s\n", filename);
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Read number of elements */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_items) != 1)
    {
        fprintf(stderr, "Error: failed to read number of elements in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }
    if (num_items != num_elements)
    {
        fprintf(stderr, "Error: invalid number of elements in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Skip elements */
    for (i = 0; i < num_elements; i++)
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to skip element %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

    /* Read number of nodes */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_items) != 1)
    {
        fprintf(stderr, "Error: failed to read number of nodes in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }
    if (num_items != num_nodes)
    {
        fprintf(stderr, "Error: invalid number of nodes in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Skip nodes */
    for (i = 0; i < num_nodes; i++)
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to skip node %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

    /* Read number of material IDs */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_material_ids) != 1)
    {
        fprintf(stderr, "Error: failed to read number of material IDs in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Allocate integer elemental fields */
    if (num_elements > previous_num_elements)
    {
        num_integer_elemental_fields = 1;

        flag = (void *)integer_elemental_fields;
        integer_elemental_fields
            = (int **)realloc(integer_elemental_fields,
                              sizeof(int *) * num_integer_elemental_fields);
        if (num_integer_elemental_fields)
        {
            if (flag == NULL) integer_elemental_fields[0] = NULL;
            integer_elemental_fields[0]
                = (int *)realloc(integer_elemental_fields[0],
                                 sizeof(int) * num_elements * num_integer_elemental_fields);
        }
        for (i = 0; i < num_integer_elemental_fields; i++)
            integer_elemental_fields[i]
                = integer_elemental_fields[0] + num_elements * i;

        casted_integer_elemental_field
            = (double *)realloc(casted_integer_elemental_field,
                                sizeof(double) * num_elements);

        flag = (void *)integer_elemental_field_labels;
        integer_elemental_field_labels
            = (char **)realloc(integer_elemental_field_labels,
                               sizeof(char *) * num_integer_elemental_fields);
        if (num_integer_elemental_fields)
        {
            if (flag == NULL) integer_elemental_field_labels[0] = NULL;
            integer_elemental_field_labels[0]
                = (char *)realloc(integer_elemental_field_labels[0],
                                  sizeof(char) * 256 * num_integer_elemental_fields);
        }
        for (i = 0; i < num_integer_elemental_fields; i++)
            integer_elemental_field_labels[i]
                = integer_elemental_field_labels[0] + 256 * i;

        previous_num_elements = num_elements;
    }

    /* Set label */
    strcpy(integer_elemental_field_labels[0], "MaterialID");

    /* Read material IDs */
    for (i = 0; i < num_material_ids; i++)
    {
        int num_material_elements;
        int material_element_id;

        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer, "%d", &num_material_elements) != 1)
        {
            fprintf(stderr, "Error: failed to read number of material ID %d elements in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < num_material_elements; j++)
        {
            if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
                || sscanf(buffer, "%d", &material_element_id) != 1)
            {
                fprintf(stderr, "Error: failed to read material ID %d element ID %d in %s\n",
                        i, j, filename);
                exit(EXIT_FAILURE);
            }

            integer_elemental_fields[0][material_element_id] = i;
        }
    }

    /* Close */
    fclose(fp);
}

/*
 * Init or reinit fields from emsh file
 */
static void initOrReinitFieldsFromEmshFile(void)
{
    FILE *fp;

    char filename[256];
    char buffer[FILE_INPUT_BUFFER_SIZE];

    void *flag;

    static int previous_num_elements = -1;

    int num_items;
    int num_material_ids;
    int i, j;

    /* Open file */
    sprintf(filename, filenames_fields[0],
            current_time_step, current_time_step, current_time_step);
    printf("Input %s\n", filename);
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Read number of elements */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_items) != 1)
    {
        fprintf(stderr, "Error: failed to read number of elements in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }
    if (num_items != num_elements)
    {
        fprintf(stderr, "Error: invalid number of elements in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Skip elements */
    for (i = 0; i < num_elements; i++)
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to skip element %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

    /* Read number of segments */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_items) != 1)
    {
        fprintf(stderr, "Error: failed to read number of segments in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }
    if (num_items != num_segments)
    {
        fprintf(stderr, "Error: invalid number of segments in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Skip segments */
    for (i = 0; i < num_segments; i++)
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to skip segment %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

    /* Read number of nodes */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_items) != 1)
    {
        fprintf(stderr, "Error: failed to read number of nodes in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }
    if (num_items != num_nodes)
    {
        fprintf(stderr, "Error: invalid number of nodes in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Skip nodes */
    for (i = 0; i < num_nodes; i++)
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to skip node %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

    /* Read number of material IDs */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_material_ids) != 1)
    {
        fprintf(stderr, "Error: failed to read number of material IDs in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Allocate integer elemental fields */
    if (num_elements > previous_num_elements)
    {
        num_integer_elemental_fields = 1;

        flag = (void *)integer_elemental_fields;
        integer_elemental_fields
            = (int **)realloc(integer_elemental_fields,
                              sizeof(int *) * num_integer_elemental_fields);
        if (num_integer_elemental_fields)
        {
            if (flag == NULL) integer_elemental_fields[0] = NULL;
            integer_elemental_fields[0]
                = (int *)realloc(integer_elemental_fields[0],
                                 sizeof(int) * num_elements * num_integer_elemental_fields);
        }
        for (i = 0; i < num_integer_elemental_fields; i++)
            integer_elemental_fields[i]
                = integer_elemental_fields[0] + num_elements * i;

        casted_integer_elemental_field
            = (double *)realloc(casted_integer_elemental_field,
                                sizeof(double) * num_elements);

        flag = (void *)integer_elemental_field_labels;
        integer_elemental_field_labels
            = (char **)realloc(integer_elemental_field_labels,
                               sizeof(char *) * num_integer_elemental_fields);
        if (num_integer_elemental_fields)
        {
            if (flag == NULL) integer_elemental_field_labels[0] = NULL;
            integer_elemental_field_labels[0]
                = (char *)realloc(integer_elemental_field_labels[0],
                                  sizeof(char) * 256 * num_integer_elemental_fields);
        }
        for (i = 0; i < num_integer_elemental_fields; i++)
            integer_elemental_field_labels[i]
                = integer_elemental_field_labels[0] + 256 * i;

        previous_num_elements = num_elements;
    }

    /* Set label */
    strcpy(integer_elemental_field_labels[0], "MaterialID");

    /* Read material IDs */
    for (i = 0; i < num_material_ids; i++)
    {
        int num_material_elements;
        int material_element_id;

        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer, "%d", &num_material_elements) != 1)
        {
            fprintf(stderr, "Error: failed to read number of material ID %d elements in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

        for (j = 0; j < num_material_elements; j++)
        {
            if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
                || sscanf(buffer, "%d", &material_element_id) != 1)
            {
                fprintf(stderr, "Error: failed to read material ID %d element ID %d in %s\n",
                        i, j, filename);
                exit(EXIT_FAILURE);
            }

            integer_elemental_fields[0][material_element_id] = i;
        }
    }

    /* Close */
    fclose(fp);
}

/*
 * Init or reinit fields from pcm file
 */
static void initOrReinitFieldsFromPcmFile(void)
{
    FILE *fp;

    char filename[256];
    char buffer[FILE_INPUT_BUFFER_SIZE];

    void *flag;

    static int previous_num_elements = -1;

    int num_items;
    int num_material_ids;
    int temp1, temp2;
    int offset;
    int i, j;

    /* Open file */
    sprintf(filename, filenames_fields[0],
            current_time_step, current_time_step, current_time_step);
    printf("Input %s\n", filename);
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Read numbers of nodes and material IDs */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d%d%d",
                  &num_items, &temp1, &num_material_ids) != 3)
    {
        fprintf(stderr, "Error: failed to read numbers of nodes and material IDs in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }
    if (num_items != num_nodes)
    {
        fprintf(stderr, "Error: invalid number of nodes in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Skip nodes */
    for (i = 0; i < num_nodes; i++)
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to skip node %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

    /* Allocate integer elemental fields */
    if (num_elements > previous_num_elements)
    {
        num_integer_elemental_fields = 1;

        flag = (void *)integer_elemental_fields;
        integer_elemental_fields
            = (int **)realloc(integer_elemental_fields,
                              sizeof(int *) * num_integer_elemental_fields);
        if (num_integer_elemental_fields)
        {
            if (flag == NULL) integer_elemental_fields[0] = NULL;
            integer_elemental_fields[0]
                = (int *)realloc(integer_elemental_fields[0],
                                 sizeof(int) * num_elements * num_integer_elemental_fields);
        }
        for (i = 0; i < num_integer_elemental_fields; i++)
            integer_elemental_fields[i]
                = integer_elemental_fields[0] + num_elements * i;

        casted_integer_elemental_field
            = (double *)realloc(casted_integer_elemental_field,
                                sizeof(double) * num_elements);

        flag = (void *)integer_elemental_field_labels;
        integer_elemental_field_labels
            = (char **)realloc(integer_elemental_field_labels,
                               sizeof(char *) * num_integer_elemental_fields);
        if (num_integer_elemental_fields)
        {
            if (flag == NULL) integer_elemental_field_labels[0] = NULL;
            integer_elemental_field_labels[0]
                = (char *)realloc(integer_elemental_field_labels[0],
                                  sizeof(char) * 256 * num_integer_elemental_fields);
        }
        for (i = 0; i < num_integer_elemental_fields; i++)
            integer_elemental_field_labels[i]
                = integer_elemental_field_labels[0] + 256 * i;

        previous_num_elements = num_elements;
    }

    /* Set label */
    strcpy(integer_elemental_field_labels[0], "MaterialID");

    offset = 0;
    for (i = 0; i < num_material_ids; i++)
    {
        /* Read number of elements */
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer, "%d%d%d",
                      &num_items, &temp1, &temp2) != 3)
        {
            fprintf(stderr, "Error: failed to read number of elements in %s\n",
                    filename);
            exit(EXIT_FAILURE);
        }

        /* Skip elements */
        for (j = 0; j < num_items; j++)
        {
            if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
            {
                fprintf(stderr, "Error: failed to skip element %d-%d in %s\n",
                        i, j, filename);
                exit(EXIT_FAILURE);
            }

            integer_elemental_fields[0][offset++] = i;
        }
    }
    if (offset != num_elements)
    {
        fprintf(stderr, "Error: invalid number of elements in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Close */
    fclose(fp);
}

/*
 * Init or reinit fields from visible elements file
 */
static void initOrReinitFieldsFromVisibleElementsFile(void)
{
    FILE *fp;

    char filename[256];
    char buffer[FILE_INPUT_BUFFER_SIZE];

    void *flag;

    static int previous_num_elements = -1;

    int num_items;
    int element_id;
    int i;

    /* Open file */
    sprintf(filename, filenames_fields[0],
            current_time_step, current_time_step, current_time_step);
    printf("Input %s\n", filename);
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    /* Allocate integer elemental fields */
    if (num_elements > previous_num_elements)
    {
        num_integer_elemental_fields = 1;

        flag = (void *)integer_elemental_fields;
        integer_elemental_fields
            = (int **)realloc(integer_elemental_fields,
                              sizeof(int *) * num_integer_elemental_fields);
        if (num_integer_elemental_fields)
        {
            if (flag == NULL) integer_elemental_fields[0] = NULL;
            integer_elemental_fields[0]
                = (int *)realloc(integer_elemental_fields[0],
                                 sizeof(int) * num_elements * num_integer_elemental_fields);
        }
        for (i = 0; i < num_integer_elemental_fields; i++)
            integer_elemental_fields[i]
                = integer_elemental_fields[0] + num_elements * i;

        casted_integer_elemental_field
            = (double *)realloc(casted_integer_elemental_field,
                                sizeof(double) * num_elements);

        flag = (void *)integer_elemental_field_labels;
        integer_elemental_field_labels
            = (char **)realloc(integer_elemental_field_labels,
                               sizeof(char *) * num_integer_elemental_fields);
        if (num_integer_elemental_fields)
        {
            if (flag == NULL) integer_elemental_field_labels[0] = NULL;
            integer_elemental_field_labels[0]
                = (char *)realloc(integer_elemental_field_labels[0],
                                  sizeof(char) * 256 * num_integer_elemental_fields);
        }
        for (i = 0; i < num_integer_elemental_fields; i++)
            integer_elemental_field_labels[i]
                = integer_elemental_field_labels[0] + 256 * i;

        previous_num_elements = num_elements;
    }

    /* Set label */
    strcpy(integer_elemental_field_labels[0], "VisibleElements");

    /* Fill integer elemental fields */
    for (i = 0; i < num_elements; i++)
        integer_elemental_fields[0][i] = 0;

    /* Read visible elements */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_items) != 1)
    {
        fprintf(stderr, "Error: failed to read number of visible elements in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < num_items; i++)
    {
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer, "%d", &element_id) != 1)
        {
            fprintf(stderr, "Error: failed to read visible element %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

        integer_elemental_fields[0][element_id] = 1;
    }

    /* Close */
    fclose(fp);
}

/*
 * Destroy fields
 */
void destroyFields(void)
{
    if (num_displacement_fields)
    {
        free(displacement_field_labels[0]);
        free(displacement_field_labels);
        free(displacement_field_offsets);
    }

    if (num_integer_elemental_fields)
        free(integer_elemental_field_labels[0]);
    free(integer_elemental_field_labels);
    free(casted_integer_elemental_field);
    if (num_integer_elemental_fields)
        free(integer_elemental_fields[0]);
    free(integer_elemental_fields);

    if (num_elemental_fields)
        free(elemental_field_labels[0]);
    free(elemental_field_labels);
    if (num_elemental_fields)
        free(elemental_fields[0]);
    free(elemental_fields);

    if (num_nodal_fields)
        free(nodal_field_labels[0]);
    free(nodal_field_labels);
    if (num_nodal_fields)
        free(nodal_fields[0]);
    free(nodal_fields);

    free(zero_vector);
}

/*
 * Set contour field
 */
void setContourField(void)
{
    int i;

    if (contour_field_id >= 0
        && contour_field_id < num_nodal_fields)
    {
        int nodal_field_id = contour_field_id;
        contour_field = nodal_fields[nodal_field_id];
        strcpy(contour_bar_label,
               nodal_field_labels[nodal_field_id]);
    }
    else if (contour_field_id >= num_nodal_fields
             && contour_field_id < num_nodal_fields + num_elemental_fields)
    {
        int elemental_field_id = contour_field_id - num_nodal_fields;
        contour_field = elemental_fields[elemental_field_id];
        strcpy(contour_bar_label,
               elemental_field_labels[elemental_field_id]);
    }
    else if (contour_field_id >= num_nodal_fields + num_elemental_fields
             && contour_field_id < num_nodal_fields + num_elemental_fields + num_integer_elemental_fields)
    {
        int integer_elemental_field_id = contour_field_id - num_nodal_fields - num_elemental_fields;
        contour_field = casted_integer_elemental_field;
        strcpy(contour_bar_label,
               integer_elemental_field_labels[integer_elemental_field_id]);

        if (integer_contour_field_mod == 0)
        {
            for (i = 0; i < num_elements; i++)
                casted_integer_elemental_field[i]
                    = (double)integer_elemental_fields[integer_elemental_field_id][i];
        }
        else
        {
            for (i = 0; i < num_elements; i++)
                casted_integer_elemental_field[i]
                    = (double)(integer_elemental_fields[integer_elemental_field_id][i] % integer_contour_field_mod);
        }
    }
    else
    {
        contour_field = NULL;
        strcpy(contour_bar_label, "");
    }
}

/*
 * Set contour field range
 */
void setContourFieldRange(void)
{
    int num;
    int i;

    if (contour_field_id >= 0
        && contour_field_id < num_nodal_fields)
        num = num_nodes;
    else if (contour_field_id >= num_nodal_fields
             && contour_field_id < num_nodal_fields + num_elemental_fields)
        num = num_elements;
    else if (contour_field_id >= num_nodal_fields + num_elemental_fields
             && contour_field_id < num_nodal_fields + num_elemental_fields + num_integer_elemental_fields)
        num = num_elements;
    else
        num = 0;

    if (contour_field)
    {
        contour_field_max = -DBL_MAX;
        contour_field_min = DBL_MAX;

        for (i = 0; i < num; i++)
        {
            if (contour_field_max < contour_field[i])
                contour_field_max = contour_field[i];
            if (contour_field_min > contour_field[i])
                contour_field_min = contour_field[i];
        }
    }
    else
    {
        contour_field_max = 1.0;
        contour_field_min = 0.0;
    }
}

/*
 * Set displacement field
 */
void setDisplacementField(void)
{
    if (displacement_field_id >= 0
        && displacement_field_id + num_dimensions <= num_nodal_fields)
    {
        int offset = displacement_field_offsets[displacement_field_id];

        displacement_fields[0] = nodal_fields[offset];
        displacement_fields[1] = nodal_fields[offset + 1];
        if (num_dimensions == 3)
            displacement_fields[2] = nodal_fields[offset + 2];
        else
            displacement_fields[2] = zero_vector;
    }
    else
    {
        displacement_fields[0] = NULL;
        displacement_fields[1] = NULL;
        displacement_fields[2] = NULL;
    }
}
