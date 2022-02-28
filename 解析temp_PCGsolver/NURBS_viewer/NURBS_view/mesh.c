/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef WITH_ADVIO
#include <Adv/AdvDocument.h>
#endif

#include "mesh.h"

#include "constant.h"
#include "face.h"
#include "edge.h"
#include "crack.h"

extern char *filename_mesh;
extern char filename_visible_elements[];
extern char filename_invisible_elements[];

extern int current_time_step;

extern int element_type;
extern int num_elements;
extern int num_element_nodes;
extern int **elements;
extern int num_visible_elements;
extern int *visible_elements;
extern int num_visible_tetrahedrons;
extern int num_visible_hexahedrons;
extern int *visible_tetrahedrons;
extern int *visible_hexahedrons;
extern int num_segments;
extern int **segments;
extern int num_nodes;
extern int num_dimensions;
extern double **nodes;

static void initOrReinitMeshFromAdvFile(void);
static void initOrReinitMeshFromMshFile(int extended_mesh_flag);
static void initOrReinitMeshFromPchFile(void);
static void initOrReinitMeshFromPcmFile(void);
static void initOrReinitVisibleElements(void);
static void reinitVisibleElementsFromVisibleAndInvisibleElementsFile(void);
static void initOrReinitVisibleTetrahedronsAndHexahedrons(void);

/*
 * Init or reinit mesh
 */
void initOrReinitMesh(void)
{
    char filename[256];

    static char previous_filename[256] = "";

    char *filename_suffix;

    /* Check filename */
    sprintf(filename, filename_mesh,
            current_time_step, current_time_step, current_time_step);
    if (strcmp(filename, previous_filename) == 0) return;
    strcpy(previous_filename, filename);

    /* Get filename suffix */
    filename_suffix = strrchr(filename_mesh, '.');
    if (filename_suffix == NULL)
        filename_suffix = filename_mesh;

    /* Init or reinit mesh */
    if (strcmp(filename_suffix, ".adv") == 0)
        initOrReinitMeshFromAdvFile();
    else if (strcmp(filename_suffix, ".msh") == 0)
        initOrReinitMeshFromMshFile(0);
    else if (strcmp(filename_suffix, ".emsh") == 0)
        initOrReinitMeshFromMshFile(1);
    else if (strcmp(filename_suffix, ".pch") == 0)
        initOrReinitMeshFromPchFile();
    else if (strcmp(filename_suffix, ".pcm") == 0)
        initOrReinitMeshFromPcmFile();
    else
        initOrReinitMeshFromMshFile(0);

    /* Init or reinit visible elements */
    initOrReinitVisibleElements();

    /* Init or reinit visible tetrahedrons and hexahedrons */
    initOrReinitVisibleTetrahedronsAndHexahedrons();

    /* Init or reinit faces */
    initOrReinitFaces();

    /* Init or reinit edges */
    initOrReinitEdges();

    /* Init or reinit crack faces */
    initOrReinitCrackFaces();

    /* Init or reinit crack front edges */
    initOrReinitCrackFrontEdges();
}

/*
 * Init or reinit mesh from adv file
 */
static void initOrReinitMeshFromAdvFile(void)
#ifdef WITH_ADVIO
{
    AdvDocFile *dfile;
    AdvDocument *doc;
    adv_off_t offset;

    char filename[256];

    void *flag;

    static int previous_num_elements = -1;
    static int previous_num_element_nodes = -1;
    static int previous_num_nodes = -1;

    char element_type_property[256];
    int i;

    /* Open file */
    sprintf(filename, filename_mesh,
            current_time_step, current_time_step, current_time_step);
    printf("Input %s\n", filename);
    dfile = adv_dio_file_open(filename, "r");
    if (dfile == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    /* Open elements */
    doc = adv_dio_open_by_property(dfile, NULL,
                                   "content_type", "Element",
                                   NULL);
    if (doc == NULL)
    {
        fprintf(stderr, "Error: failed to open document %s in %s\n",
                "Element", filename);
        exit(EXIT_FAILURE);
    }

    /* Get properties */
    strcpy(element_type_property, adv_dio_get_property(doc, "element_type"));
    adv_dio_get_property_int32(doc, "num_items", &num_elements);
    adv_dio_get_property_int32(doc, "num_nodes_per_element", &num_element_nodes);

    printf("Element type: %s\n", element_type_property);
    printf("Number of elements: %d\n", num_elements);
    printf("Number of element nodes: %d\n", num_element_nodes);

    /* Allocate elements */
    if (num_elements > previous_num_elements
        || num_element_nodes > previous_num_element_nodes)
    {
        flag = (void *)elements;
        elements
            = (int **)realloc(elements,
                              sizeof(int *) * num_elements);
        if (num_elements)
        {
            if (flag == NULL) elements[0] = NULL;
            elements[0]
                = (int *)realloc(elements[0],
                                 sizeof(int) * num_element_nodes * num_elements);
        }
        for (i = 0; i < num_elements; i++)
            elements[i]
                = elements[0] + num_element_nodes * i;

        previous_num_elements = num_elements;
        previous_num_element_nodes = num_element_nodes;
    }

    /* Read elements */
    adv_dio_read_int32v(doc, 0,
                        num_element_nodes * num_elements,
                        &elements[0][0]);

    /* Close elements */
    adv_dio_close(doc);

    /* Open nodes */
    doc = adv_dio_open_by_property(dfile, NULL,
                                   "content_type", "Node",
                                   NULL);
    if (doc == NULL)
    {
        fprintf(stderr, "Error: failed to open document %s in %s\n",
                "Node", filename);
        exit(EXIT_FAILURE);
    }

    /* Get properties */
    adv_dio_get_property_int32(doc, "num_items", &num_nodes);
    adv_dio_get_property_int32(doc, "dimension", &num_dimensions);
    if (num_dimensions != 3 && num_dimensions != 2)
    {
        fprintf(stderr, "Error: unsupported number of dimensions (%d)\n",
                num_dimensions);
        exit(EXIT_FAILURE);
    }

    printf("Number of nodes: %d\n", num_nodes);
    printf("Number of dimensions: %d\n", num_dimensions);

    /* Allocate nodes */
    if (num_nodes > previous_num_nodes)
    {
        flag = (void *)nodes;
        nodes
            = (double **)realloc(nodes,
                                 sizeof(double *) * num_nodes);
        if (num_nodes)
        {
            if (flag == NULL) nodes[0] = NULL;
            nodes[0]
                = (double *)realloc(nodes[0],
                                    sizeof(double) * 3 * num_nodes);
        }
        for (i = 0; i < num_nodes; i++)
            nodes[i]
                = nodes[0] + 3 * i;

        previous_num_nodes = num_nodes;
    }

    /* Read nodes */
    offset = 0;
    for (i = 0; i < num_nodes; i++)
        offset += adv_dio_read_float64v(doc, offset,
                                        num_dimensions,
                                        nodes[i]);

    /* Close nodes */
    adv_dio_close(doc);

    /* Close */
    adv_dio_file_close(dfile);

    /* Zero-fill nodal coordinate Z if 2D */
    if (num_dimensions == 2)
        for (i = 0; i < num_nodes; i++)
            nodes[i][2] = 0.0;

    /* Get element type */
    if (strcmp(element_type_property, "3DLinearTetrahedron") == 0)
        element_type = LINEAR_TETRAHEDRON_ELEMENT_TYPE;
    else if (strcmp(element_type_property, "3DQuadraticTetrahedron") == 0)
        element_type = QUADRATIC_TETRAHEDRON_ELEMENT_TYPE;
    else if (strcmp(element_type_property, "3DLinearHexahedron") == 0)
        element_type = LINEAR_HEXAHEDRON_ELEMENT_TYPE;
    else if (strcmp(element_type_property, "3DQuadraticHexahedron") == 0)
        element_type = QUADRATIC_HEXAHEDRON_ELEMENT_TYPE;
    else if (strcmp(element_type_property, "3DQuadraticHexahedron27") == 0)
        element_type = QUADRATIC_HEXAHEDRON27_ELEMENT_TYPE;
    else if (strcmp(element_type_property, "3DLinearTriangle") == 0
             || strcmp(element_type_property, "2DLinearTriangle") == 0)
        element_type = LINEAR_TRIANGLE_ELEMENT_TYPE;
    else if (strcmp(element_type_property, "3DQuadraticTriangle") == 0
             || strcmp(element_type_property, "2DQuadraticTriangle") == 0)
        element_type = QUADRATIC_TRIANGLE_ELEMENT_TYPE;
    else if (strcmp(element_type_property, "3DLinearQuadrangle") == 0
             || strcmp(element_type_property, "2DLinearQuadrangle") == 0)
        element_type = LINEAR_QUADRANGLE_ELEMENT_TYPE;
    else if (strcmp(element_type_property, "3DQuadraticQuadrangle") == 0
             || strcmp(element_type_property, "2DQuadraticQuadrangle") == 0)
        element_type = QUADRATIC_QUADRANGLE_ELEMENT_TYPE;
    else if (strcmp(element_type_property, "3DQuadraticQuadrangle2") == 0
             || strcmp(element_type_property, "2DQuadraticQuadrangle2") == 0)
        element_type = QUADRATIC_QUADRANGLE9_ELEMENT_TYPE;
    else
    {
        fprintf(stderr, "Error: unsupported element type %s\n",
                element_type_property);
        exit(EXIT_FAILURE);
    }
}
#else
{
    fprintf(stderr, "Error: adv file not supported in this build\n");
    exit(EXIT_FAILURE);
}
#endif

/*
 * Init or reinit mesh from msh file
 */
static void initOrReinitMeshFromMshFile(int extended_mesh_flag)
{
    FILE *fp;
    fpos_t position;

    char filename[256];
    char buffer[FILE_INPUT_BUFFER_SIZE];

    void *flag;

    static int previous_num_elements = -1;
    static int previous_num_element_nodes = -1;
    static int previous_num_segments = -1;
    static int previous_num_nodes = -1;

    int element_nodes[27];
    double nodal_coords[3];
    int i;

    /* Open file */
    sprintf(filename, filename_mesh,
            current_time_step, current_time_step, current_time_step);
    printf("Input %s\n", filename);
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    /* Read number of elements */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_elements) != 1)
    {
        fprintf(stderr, "Error: failed to read number of elements in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    printf("Number of elements: %d\n", num_elements);

    /* Read 1st element */
    if (num_elements)
    {
        fgetpos(fp, &position);
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to read element 0 in %s\n", filename);
            exit(EXIT_FAILURE);
        }
        fsetpos(fp, &position);
        num_element_nodes = sscanf(buffer,
                                   "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
                                   &element_nodes[0],
                                   &element_nodes[1],
                                   &element_nodes[2],
                                   &element_nodes[3],
                                   &element_nodes[4],
                                   &element_nodes[5],
                                   &element_nodes[6],
                                   &element_nodes[7],
                                   &element_nodes[8],
                                   &element_nodes[9],
                                   &element_nodes[10],
                                   &element_nodes[11],
                                   &element_nodes[12],
                                   &element_nodes[13],
                                   &element_nodes[14],
                                   &element_nodes[15],
                                   &element_nodes[16],
                                   &element_nodes[17],
                                   &element_nodes[18],
                                   &element_nodes[19],
                                   &element_nodes[20],
                                   &element_nodes[21],
                                   &element_nodes[22],
                                   &element_nodes[23],
                                   &element_nodes[24],
                                   &element_nodes[25],
                                   &element_nodes[26]);
    }
    else
    {
        num_element_nodes = 3;
    }

    printf("Number of element nodes: %d\n", num_element_nodes);

    /* Allocate elements */
    if (num_elements > previous_num_elements
        || num_element_nodes > previous_num_element_nodes)
    {
        flag = (void *)elements;
        elements
            = (int **)realloc(elements,
                              sizeof(int *) * num_elements);
        if (num_elements)
        {
            if (flag == NULL) elements[0] = NULL;
            elements[0]
                = (int *)realloc(elements[0],
                                 sizeof(int) * num_element_nodes * num_elements);
        }
        for (i = 0; i < num_elements; i++)
            elements[i]
                = elements[0] + num_element_nodes * i;

        previous_num_elements = num_elements;
        previous_num_element_nodes = num_element_nodes;
    }

    /* Read elements */
    for (i = 0; i < num_elements; i++)
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer,
                      "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
                      &elements[i][0],
                      &elements[i][1],
                      &elements[i][2],
                      &elements[i][3],
                      &elements[i][4],
                      &elements[i][5],
                      &elements[i][6],
                      &elements[i][7],
                      &elements[i][8],
                      &elements[i][9],
                      &elements[i][10],
                      &elements[i][11],
                      &elements[i][12],
                      &elements[i][13],
                      &elements[i][14],
                      &elements[i][15],
                      &elements[i][16],
                      &elements[i][17],
                      &elements[i][18],
                      &elements[i][19],
                      &elements[i][20],
                      &elements[i][21],
                      &elements[i][22],
                      &elements[i][23],
                      &elements[i][24],
                      &elements[i][25],
                      &elements[i][26]) != num_element_nodes)
        {
            fprintf(stderr, "Error: failed to read element %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

    if (extended_mesh_flag)
    {
        /* Read number of segments */
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer, "%d", &num_segments) != 1)
        {
            fprintf(stderr, "Error: failed to read number of segments in %s\n",
                    filename);
            exit(EXIT_FAILURE);
        }

        printf("Number of segments: %d\n", num_segments);

        /* Allocate segments */
        if (num_segments > previous_num_segments)
        {
            flag = (void *)segments;
            segments
                = (int **)realloc(segments,
                                  sizeof(int *) * num_segments);
            if (num_segments)
            {
                if (flag == NULL) segments[0] = NULL;
                segments[0]
                    = (int *)realloc(segments[0],
                                     sizeof(int) * 2 * num_segments);
            }
            for (i = 0; i < num_segments; i++)
                segments[i]
                    = segments[0] + 2 * i;

            previous_num_segments = num_segments;
        }

        /* Read segments */
        for (i = 0; i < num_segments; i++)
            if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
                || sscanf(buffer,
                          "%d%d",
                          &segments[i][0],
                          &segments[i][1]) != 2)
            {
                fprintf(stderr, "Error: failed to read segment %d in %s\n",
                        i, filename);
                exit(EXIT_FAILURE);
            }
    }

    /* Read number of nodes */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_nodes) != 1)
    {
        fprintf(stderr, "Error: failed to read number of nodes in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    printf("Number of nodes: %d\n", num_nodes);

    /* Read 1st node */
    if (num_nodes)
    {
        fgetpos(fp, &position);
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to read node 0 in %s\n", filename);
            exit(EXIT_FAILURE);
        }
        fsetpos(fp, &position);
        num_dimensions = sscanf(buffer,
                                "%lf%lf%lf",
                                &nodal_coords[0],
                                &nodal_coords[1],
                                &nodal_coords[2]);
        if (num_dimensions != 3 && num_dimensions != 2)
        {
            fprintf(stderr, "Error: unsupported number of dimensions (%d)\n",
                    num_dimensions);
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        num_dimensions = 3;
    }

    printf("Number of dimensions: %d\n", num_dimensions);

    /* Allocate nodes */
    if (num_nodes > previous_num_nodes)
    {
        flag = (void *)nodes;
        nodes
            = (double **)realloc(nodes,
                                 sizeof(double *) * num_nodes);
        if (num_nodes)
        {
            if (flag == NULL) nodes[0] = NULL;
            nodes[0]
                = (double *)realloc(nodes[0],
                                    sizeof(double) * 3 * num_nodes);
        }
        for (i = 0; i < num_nodes; i++)
            nodes[i]
                = nodes[0] + 3 * i;

        previous_num_nodes = num_nodes;
    }

    /* Read nodes */
    for (i = 0; i < num_nodes; i++)
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer,
                      "%lf%lf%lf",
                      &nodes[i][0],
                      &nodes[i][1],
                      &nodes[i][2]) != num_dimensions)
        {
            fprintf(stderr, "Error: failed to read node %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

    /* Close */
    fclose(fp);

    /* Zero-fill nodal coordinate Z if 2D */
    if (num_dimensions == 2)
        for (i = 0; i < num_nodes; i++)
            nodes[i][2] = 0.0;

    /* Get element type */
    if (num_dimensions == 3)
        switch (num_element_nodes)
        {
        case 4:
            element_type = LINEAR_TETRAHEDRON_ELEMENT_TYPE;
            break;
        case 10:
            element_type = QUADRATIC_TETRAHEDRON_ELEMENT_TYPE;
            break;
        case 8:
            element_type = LINEAR_HEXAHEDRON_ELEMENT_TYPE;
            break;
        case 20:
            element_type = QUADRATIC_HEXAHEDRON_ELEMENT_TYPE;
            break;
        case 27:
            element_type = QUADRATIC_HEXAHEDRON27_ELEMENT_TYPE;
            break;
        case 3:
            element_type = LINEAR_TRIANGLE_ELEMENT_TYPE;
            break;
        case 6:
            element_type = QUADRATIC_TRIANGLE_ELEMENT_TYPE;
            break;
            /* case 4:
             *     element_type = LINEAR_QUADRANGLE_ELEMENT_TYPE;
             *     break; */
            /* case 8:
             *     element_type = QUADRATIC_QUADRANGLE_ELEMENT_TYPE;
             *     break; */
            case 9:
                element_type = QUADRATIC_QUADRANGLE9_ELEMENT_TYPE;
                break;
        default:
            fprintf(stderr, "Error: unsupported element type (%d element nodes)\n",
                    num_element_nodes);
            exit(EXIT_FAILURE);
        }
    else
        switch (num_element_nodes)
        {
        case 3:
            element_type = LINEAR_TRIANGLE_ELEMENT_TYPE;
            break;
        case 6:
            element_type = QUADRATIC_TRIANGLE_ELEMENT_TYPE;
            break;
        case 4:
            element_type = LINEAR_QUADRANGLE_ELEMENT_TYPE;
            break;
        case 8:
            element_type = QUADRATIC_QUADRANGLE_ELEMENT_TYPE;
            break;
        case 9:
            element_type = QUADRATIC_QUADRANGLE9_ELEMENT_TYPE;
            break;
        default:
            fprintf(stderr, "Error: unsupported element type (%d element nodes)\n",
                    num_element_nodes);
            exit(EXIT_FAILURE);
        }
}

/*
 * Init or reinit mesh from pch file
 */
static void initOrReinitMeshFromPchFile(void)
{
    FILE *fp;
    fpos_t position;

    char filename[256];
    char buffer[FILE_INPUT_BUFFER_SIZE];

    void *flag;

    static int previous_num_nodes = -1;
    static int previous_num_elements = -1;
    static int previous_num_element_nodes = -1;

    double nodal_coords[3];
    int element_nodes[8];
    int i;

    /* Open file */
    sprintf(filename, filename_mesh,
            current_time_step, current_time_step, current_time_step);
    printf("Input %s\n", filename);
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    /* Read number of nodes */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_nodes) != 1)
    {
        fprintf(stderr, "Error: failed to read number of nodes in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    printf("Number of nodes: %d\n", num_nodes);

    /* Read 1st node */
    if (num_nodes)
    {
        fgetpos(fp, &position);
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to read node 0 in %s\n", filename);
            exit(EXIT_FAILURE);
        }
        fsetpos(fp, &position);
        num_dimensions = sscanf(buffer,
                                "%lf%lf%lf",
                                &nodal_coords[0],
                                &nodal_coords[1],
                                &nodal_coords[2]);
        if (num_dimensions != 3 && num_dimensions != 2)
        {
            fprintf(stderr, "Error: unsupported number of dimensions (%d)\n",
                    num_dimensions);
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        num_dimensions = 3;
    }

    printf("Number of dimensions: %d\n", num_dimensions);

    /* Allocate nodes */
    if (num_nodes > previous_num_nodes)
    {
        flag = (void *)nodes;
        nodes
            = (double **)realloc(nodes,
                                 sizeof(double *) * num_nodes);
        if (num_nodes)
        {
            if (flag == NULL) nodes[0] = NULL;
            nodes[0]
                = (double *)realloc(nodes[0],
                                    sizeof(double) * 3 * num_nodes);
        }
        for (i = 0; i < num_nodes; i++)
            nodes[i]
                = nodes[0] + 3 * i;

        previous_num_nodes = num_nodes;
    }

    /* Read nodes */
    for (i = 0; i < num_nodes; i++)
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer,
                      "%lf%lf%lf",
                      &nodes[i][0],
                      &nodes[i][1],
                      &nodes[i][2]) != num_dimensions)
        {
            fprintf(stderr, "Error: failed to read node %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

    /* Read number of elements */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d", &num_elements) != 1)
        num_elements = 0;

    printf("Number of elements: %d\n", num_elements);

    /* Read 1st element */
    if (num_elements)
    {
        fgetpos(fp, &position);
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to read element 0 in %s\n", filename);
            exit(EXIT_FAILURE);
        }
        fsetpos(fp, &position);
        num_element_nodes = sscanf(buffer,
                                   "%d%d%d%d%d%d%d%d",
                                   &element_nodes[0],
                                   &element_nodes[1],
                                   &element_nodes[2],
                                   &element_nodes[3],
                                   &element_nodes[4],
                                   &element_nodes[5],
                                   &element_nodes[6],
                                   &element_nodes[7]);
    }
    else
    {
        num_element_nodes = 3;
    }

    printf("Number of element nodes: %d\n", num_element_nodes);

    /* Allocate elements */
    if (num_elements > previous_num_elements
        || num_element_nodes > previous_num_element_nodes)
    {
        flag = (void *)elements;
        elements
            = (int **)realloc(elements,
                              sizeof(int *) * num_elements);
        if (num_elements)
        {
            if (flag == NULL) elements[0] = NULL;
            elements[0]
                = (int *)realloc(elements[0],
                                 sizeof(int) * num_element_nodes * num_elements);
        }
        for (i = 0; i < num_elements; i++)
            elements[i]
                = elements[0] + num_element_nodes * i;

        previous_num_elements = num_elements;
        previous_num_element_nodes = num_element_nodes;
    }

    /* Read elements */
    for (i = 0; i < num_elements; i++)
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer,
                      "%d%d%d%d%d%d%d%d",
                      &elements[i][0],
                      &elements[i][1],
                      &elements[i][2],
                      &elements[i][3],
                      &elements[i][4],
                      &elements[i][5],
                      &elements[i][6],
                      &elements[i][7]) != num_element_nodes)
        {
            fprintf(stderr, "Error: failed to read element %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

    /* Close */
    fclose(fp);

    /* Zero-fill nodal coordinate Z if 2D */
    if (num_dimensions == 2)
        for (i = 0; i < num_nodes; i++)
            nodes[i][2] = 0.0;

    /* Get element type */
    switch (num_element_nodes)
    {
    case 3:
        element_type = LINEAR_TRIANGLE_ELEMENT_TYPE;
        break;
    case 6:
        element_type = QUADRATIC_TRIANGLE_ELEMENT_TYPE;
        break;
    case 4:
        element_type = LINEAR_QUADRANGLE_ELEMENT_TYPE;
        break;
    case 8:
        element_type = QUADRATIC_QUADRANGLE_ELEMENT_TYPE;
        break;
    case 9:
        element_type = QUADRATIC_QUADRANGLE9_ELEMENT_TYPE;
        break;
    default:
        fprintf(stderr, "Error: unsupported element type (%d element nodes)\n",
                num_element_nodes);
        exit(EXIT_FAILURE);
    }
}

/*
 * Init or reinit mesh from pcm file
 */
static void initOrReinitMeshFromPcmFile(void)
{
    FILE *fp;
    fpos_t position;

    char filename[256];
    char buffer[FILE_INPUT_BUFFER_SIZE];

    void *flag;

    static int previous_num_nodes = -1;

    double nodal_coords[3];
    int element_nodes[8];
    int num_material_ids;
    int num_material_elements;
    int temp1, temp2;
    int i, j;

    /* Open file */
    sprintf(filename, filename_mesh,
            current_time_step, current_time_step, current_time_step);
    printf("Input %s\n", filename);
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    /* Read numbers of nodes and material IDs */
    if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
        || sscanf(buffer, "%d%d%d",
                  &num_nodes, &temp1, &num_material_ids) != 3)
    {
        fprintf(stderr, "Error: failed to read numbers of nodes and material IDs in %s\n",
                filename);
        exit(EXIT_FAILURE);
    }

    printf("Number of nodes: %d\n", num_nodes);

    /* Read 1st node */
    if (num_nodes)
    {
        fgetpos(fp, &position);
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
        {
            fprintf(stderr, "Error: failed to read node 0 in %s\n", filename);
            exit(EXIT_FAILURE);
        }
        fsetpos(fp, &position);
        num_dimensions = sscanf(buffer,
                                "%lf%lf%lf",
                                &nodal_coords[0],
                                &nodal_coords[1],
                                &nodal_coords[2]);
        if (num_dimensions != 3 && num_dimensions != 2)
        {
            fprintf(stderr, "Error: unsupported number of dimensions (%d)\n",
                    num_dimensions);
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        num_dimensions = 3;
    }

    printf("Number of dimensions: %d\n", num_dimensions);

    /* Allocate nodes */
    if (num_nodes > previous_num_nodes)
    {
        flag = (void *)nodes;
        nodes
            = (double **)realloc(nodes,
                                 sizeof(double *) * num_nodes);
        if (num_nodes)
        {
            if (flag == NULL) nodes[0] = NULL;
            nodes[0]
                = (double *)realloc(nodes[0],
                                    sizeof(double) * 3 * num_nodes);
        }
        for (i = 0; i < num_nodes; i++)
            nodes[i]
                = nodes[0] + 3 * i;

        previous_num_nodes = num_nodes;
    }

    /* Read nodes */
    for (i = 0; i < num_nodes; i++)
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer,
                      "%lf%lf%lf",
                      &nodes[i][0],
                      &nodes[i][1],
                      &nodes[i][2]) != num_dimensions)
        {
            fprintf(stderr, "Error: failed to read node %d in %s\n",
                    i, filename);
            exit(EXIT_FAILURE);
        }

    for (i = 0; i < num_material_ids; i++)
    {
        /* Read number of elements */
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer, "%d%d%d",
                      &num_material_elements, &temp1, &temp2) != 3)
        {
            fprintf(stderr, "Error: failed to read number of elements in %s\n",
                    filename);
            exit(EXIT_FAILURE);
        }
        if (i == 0) num_elements = 0;
        num_elements += num_material_elements;

        /* Read 1st element */
        if (i == 0)
        {
            if (num_material_elements)
            {
                fgetpos(fp, &position);
                if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL)
                {
                    fprintf(stderr, "Error: failed to read element 0 in %s\n",
                            filename);
                    exit(EXIT_FAILURE);
                }
                fsetpos(fp, &position);
                num_element_nodes = sscanf(buffer,
                                           "%d%d%d%d%d%d%d%d",
                                           &element_nodes[0],
                                           &element_nodes[1],
                                           &element_nodes[2],
                                           &element_nodes[3],
                                           &element_nodes[4],
                                           &element_nodes[5],
                                           &element_nodes[6],
                                           &element_nodes[7]);
            }
            else
            {
                num_element_nodes = 3;
            }
        }

        /* Allocate elements */
        flag = (void *)elements;
        elements
            = (int **)realloc(elements,
                              sizeof(int *) * num_elements);
        if (num_elements)
        {
            if (flag == NULL) elements[0] = NULL;
            elements[0]
                = (int *)realloc(elements[0],
                                 sizeof(int) * num_element_nodes * num_elements);
        }
        for (j = 0; j < num_elements; j++)
            elements[j]
                = elements[0] + num_element_nodes * j;

        /* Read elements */
        for (j = num_elements - num_material_elements;
             j < num_elements;
             j++)
            if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
                || sscanf(buffer,
                          "%d%d%d%d%d%d%d%d",
                          &elements[j][0],
                          &elements[j][1],
                          &elements[j][2],
                          &elements[j][3],
                          &elements[j][4],
                          &elements[j][5],
                          &elements[j][6],
                          &elements[j][7]) != num_element_nodes)
            {
                fprintf(stderr, "Error: failed to read element %d-%d in %s\n",
                        i, j - num_material_elements, filename);
                exit(EXIT_FAILURE);
            }
    }

    printf("Number of elements: %d\n", num_elements);
    printf("Number of element nodes: %d\n", num_element_nodes);

    /* Close */
    fclose(fp);

    /* Zero-fill nodal coordinate Z if 2D */
    if (num_dimensions == 2)
        for (i = 0; i < num_nodes; i++)
            nodes[i][2] = 0.0;

    /* Get element type */
    switch (num_element_nodes)
    {
    case 3:
        element_type = LINEAR_TRIANGLE_ELEMENT_TYPE;
        break;
    case 6:
        element_type = QUADRATIC_TRIANGLE_ELEMENT_TYPE;
        break;
    case 4:
        element_type = LINEAR_QUADRANGLE_ELEMENT_TYPE;
        break;
    case 8:
        element_type = QUADRATIC_QUADRANGLE_ELEMENT_TYPE;
        break;
    case 9:
        element_type = QUADRATIC_QUADRANGLE9_ELEMENT_TYPE;
        break;
    default:
        fprintf(stderr, "Error: unsupported element type (%d element nodes)\n",
                num_element_nodes);
        exit(EXIT_FAILURE);
    }
}

/*
 * Init or reinit visible elements
 */
static void initOrReinitVisibleElements(void)
{
    static int previous_num_visible_elements = -1;

    int i;

    /* Get number of visible elements */
    num_visible_elements = num_elements;

    /* Allocate visible elements */
    if (num_visible_elements > previous_num_visible_elements)
    {
        visible_elements
            = (int *)realloc(visible_elements,
                             sizeof(int) * num_visible_elements);

        previous_num_visible_elements = num_visible_elements;
    }

    /* Fill visible elements */
    for (i = 0; i < num_visible_elements; i++)
        visible_elements[i] = i;
}

/*
 * Destroy mesh
 */
void destroyMesh(void)
{
    destroyEdges();

    destroyCrackFaces();

    destroyFaces();

    if (num_segments)
        free(segments[0]);
    free(segments);

    free(visible_elements);

    if (num_nodes)
        free(nodes[0]);
    free(nodes);

    if (num_elements)
        free(elements[0]);
    free(elements);
}

/*
 * Reinit mesh from visible and invisible elements
 */
void reinitMeshFromVisibleAndInvisibleElements(void)
{
    char filename1[256];
    char filename2[256];

    static char previous_filename1[256] = "";
    static char previous_filename2[256] = "";

    if (strlen(filename_visible_elements) == 0
        && strlen(filename_invisible_elements) == 0)
        return;

    /* Check filenames */
    sprintf(filename1, filename_visible_elements,
            current_time_step, current_time_step, current_time_step);
    sprintf(filename2, filename_invisible_elements,
            current_time_step, current_time_step, current_time_step);
    if (strcmp(filename1, previous_filename1) == 0
        && strcmp(filename2, previous_filename2) == 0)
        return;
    strcpy(previous_filename1, filename1);
    strcpy(previous_filename2, filename2);

    /* Reinit visible elements */
    reinitVisibleElementsFromVisibleAndInvisibleElementsFile();

    /* Reinit visible tetrahedrons and hexahedrons */
    initOrReinitVisibleTetrahedronsAndHexahedrons();

    /* Reinit faces */
    reinitFacesFromVisibleElements();
}

/*
 * Reinit visible elements from visible and invisible elements file
 */
static void reinitVisibleElementsFromVisibleAndInvisibleElementsFile(void)
{
    FILE *fp;

    char filename[256];
    char buffer[FILE_INPUT_BUFFER_SIZE];

    int num_items;
    int element_id;
    int offset;
    int i;

    if (strlen(filename_visible_elements))
    {
        /* Fill visible elements */
        num_visible_elements = num_elements;
        for (i = 0; i < num_visible_elements; i++)
            visible_elements[i] = -1;

        /* Open file */
        sprintf(filename, filename_visible_elements,
                current_time_step, current_time_step, current_time_step);
        printf("Input %s\n", filename);
        fp = fopen(filename, "r");
        if (fp == NULL)
        {
            fprintf(stderr, "Error: failed to open file %s\n", filename);
            exit(EXIT_FAILURE);
        }

        /* Read number of visible elements */
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer, "%d", &num_items) != 1)
        {
            fprintf(stderr, "Error: failed to read number of visible elements in %s\n",
                    filename);
            exit(EXIT_FAILURE);
        }

        /* Read visible elements */
        for (i = 0; i < num_items; i++)
        {
            if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
                || sscanf(buffer, "%d", &element_id) != 1)
            {
                fprintf(stderr, "Error: failed to read visible element %d in %s\n",
                        i, filename);
                exit(EXIT_FAILURE);
            }

            visible_elements[element_id] = element_id;
        }

        /* Close */
        fclose(fp);
    }
    else
    {
        /* Fill visible elements */
        num_visible_elements = num_elements;
        for (i = 0; i < num_visible_elements; i++)
            visible_elements[i] = i;
    }

    if (strlen(filename_invisible_elements))
    {
        /* Open file */
        sprintf(filename, filename_invisible_elements,
                current_time_step, current_time_step, current_time_step);
        printf("Input %s\n", filename);
        fp = fopen(filename, "r");
        if (fp == NULL)
        {
            fprintf(stderr, "Error: failed to open file %s\n", filename);
            exit(EXIT_FAILURE);
        }

        /* Read number of invisible elements */
        if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
            || sscanf(buffer, "%d", &num_items) != 1)
        {
            fprintf(stderr, "Error: failed to read number of invisible elements in %s\n",
                    filename);
            exit(EXIT_FAILURE);
        }

        /* Read invisible elements */
        for (i = 0; i < num_items; i++)
        {
            if (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) == NULL
                || sscanf(buffer, "%d", &element_id) != 1)
            {
                fprintf(stderr, "Error: failed to read invisible element %d in %s\n",
                        i, filename);
                exit(EXIT_FAILURE);
            }

            visible_elements[element_id] = -1;
        }

        /* Close */
        fclose(fp);
    }

    /* Update visible elements */
    offset = 0;
    for (i = 0; i < num_visible_elements; i++)
        if (visible_elements[i] != -1)
            visible_elements[offset++] = visible_elements[i];
    num_visible_elements = offset;
}

/*
 * Init or reinit visible tetrahedrons and hexahedrons
 */
static void initOrReinitVisibleTetrahedronsAndHexahedrons(void)
{
    switch (element_type)
    {
    case LINEAR_TETRAHEDRON_ELEMENT_TYPE:
    case QUADRATIC_TETRAHEDRON_ELEMENT_TYPE:
        num_visible_tetrahedrons = num_visible_elements;
        num_visible_hexahedrons = 0;
        visible_tetrahedrons = visible_elements;
        visible_hexahedrons = NULL;
        break;
    case LINEAR_HEXAHEDRON_ELEMENT_TYPE:
    case QUADRATIC_HEXAHEDRON_ELEMENT_TYPE:
    case QUADRATIC_HEXAHEDRON27_ELEMENT_TYPE:
        num_visible_tetrahedrons = 0;
        num_visible_hexahedrons = num_visible_elements;
        visible_tetrahedrons = NULL;
        visible_hexahedrons = visible_elements;
        break;
    default:
        num_visible_tetrahedrons = 0;
        num_visible_hexahedrons = 0;
        visible_tetrahedrons = NULL;
        visible_hexahedrons = NULL;
    }
}
