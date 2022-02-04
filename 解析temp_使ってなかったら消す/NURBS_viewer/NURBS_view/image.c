#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <png.h>
#include <autogl.h>

#include "image.h"

extern char filename_view_image[];

static void writeViewImageAsPpm(FILE *fp,
                                unsigned char *image,
                                int width, int height);
static void writeViewImageAsPng(FILE *fp,
                                unsigned char *image,
                                int width, int height);

/*
 * Save view image
 */
void saveViewImage(void)
{
    FILE *fp;
    char *filename_suffix;
    unsigned char *image;
    int width, height;

    /* Allocate and get image */
    AutoGL_GetViewRangeDc(&width, &height);
    image = (unsigned char *)malloc(sizeof(unsigned char) * height * width * 3);
    AutoGL_GetViewImage(image);

    /* Get filename suffix */
    filename_suffix = strrchr(filename_view_image, '.');
    if (filename_suffix == NULL)
        filename_suffix = filename_view_image;

    /* Open file */
    fp = fopen(filename_view_image, "wb");
    if (fp == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename_view_image);
        exit(EXIT_FAILURE);
    }

    /* Write */
    if (strcmp(filename_suffix, ".ppm") == 0)
        writeViewImageAsPpm(fp, image, width, height);
    else
        writeViewImageAsPng(fp, image, width, height);

    /* Close file */
    fclose(fp);
    printf("Output %s\n", filename_view_image);

    /* Destroy image */
    free(image);
}

/*
 * Write view image as PPM
 */
static void writeViewImageAsPpm(FILE *fp,
                                unsigned char *image,
                                int width, int height)
{
    int i;

    fprintf(fp, "P6\n%d %d\n255\n", width, height);

    for (i = height - 1; i >= 0; i--)
        fwrite(&image[width * 3 * i],
               sizeof(unsigned char) * width * 3,
               sizeof(unsigned char),
               fp);
}

/*
 * Write view image as PNG
 */
static void writeViewImageAsPng(FILE *fp,
                                unsigned char *image,
                                int width, int height)
{
    png_structp write;
    png_infop info;
    png_bytep row;
    png_bytepp rows;
    int i, j;

    /* Create write structure */
    write = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (write == NULL)
    {
        fprintf(stderr, "Error: failed to create PNG write structure\n");
        exit(EXIT_FAILURE);
    }

    /* Create info structure */
    info = png_create_info_struct(write);
    if (info == NULL)
    {
        fprintf(stderr, "Error: failed to create PNG info structure\n");
        exit(EXIT_FAILURE);
    }

    /* Allocate rows */
    row  = (png_bytep)malloc(sizeof(png_byte) * height * width * 4);
    rows = (png_bytepp)malloc(sizeof(png_bytep) * height);
    for (i = 0; i < height; i++)
        rows[i] = &row[width * 4 * (height - 1 - i)];

    /* Generate rows */
    for (i = 0; i < height; i++)
        for (j = 0; j < width; j++)
        {
            png_byte red   = image[width * 3 * i + 3 * j    ];
            png_byte green = image[width * 3 * i + 3 * j + 1];
            png_byte blue  = image[width * 3 * i + 3 * j + 2];
            png_byte alpha;

            if (red == 255 && green == 255 && blue == 255)
                alpha = 0;
            else
                alpha = 255;

            row[width * 4 * i + 4 * j    ] = red;
            row[width * 4 * i + 4 * j + 1] = green;
            row[width * 4 * i + 4 * j + 2] = blue;
            row[width * 4 * i + 4 * j + 3] = alpha;
        }

    /* Write rows */
    png_init_io(write, fp);
    png_set_IHDR(write, info,
                 width, height,
                 8,
                 PNG_COLOR_TYPE_RGBA,
                 PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT);
    png_write_info(write, info);
    png_write_image(write, rows);
    png_write_end(write, info);

    /* Destroy */
    png_destroy_write_struct(&write, &info);
    free(rows);
    free(row);
}
