/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "option.h"

#include "constant.h"

extern char *filename_mesh;
extern int num_filenames_fields;
extern char **filenames_fields;
extern char filename_visible_elements[];
extern char filename_invisible_elements[];
extern char filename_points[];
extern char filename_points_fields[];
extern char filename_restart[];
extern int batch_mode_flag;
extern char filename_view_image[];
extern int window_size;
extern int division_number;

extern int first_time_step;
extern int time_step_increment;

extern int auto_view_position_flag;
extern double view_size;
extern double view_center_x;
extern double view_center_y;
extern double view_center_z;
extern double view_direction_x;
extern double view_direction_y;
extern double view_direction_z;
extern double view_up_vector_x;
extern double view_up_vector_y;
extern double view_up_vector_z;
extern double section_a;
extern double section_b;
extern double section_c;
extern double section_d;

extern int plot_type;
extern int mesh_flag;
extern int edges_flag;
extern int crack_mesh_flag;
extern int crack_front_edges_flag;
extern int segments_flag;
extern int deformed_mesh_flag;
extern int deformed_edges_flag;
extern int deformed_crack_mesh_flag;
extern int deformed_crack_front_edges_flag;
extern int deformed_segments_flag;
extern double deformation_scale_factor;
extern int points_type;
extern char title[];
extern int coord_axes_flag;
extern int coord_axes_at_origin_flag;
extern int other_mode;
extern int marked_node;
extern int marked_element;
extern int marked_point;
extern int scale_bar_flag;
extern char scale_bar_format[];
extern double scale_bar_scale_factor;
extern int section_flag;
extern int section_type;
extern char volume_color[];
extern char mesh_color[];
extern char edge_color[];
extern char crack_mesh_color[];
extern char crack_front_edge_color[];
extern char segment_color[];
extern char point_color[];
extern int mesh_line_width;
extern int edge_line_width;
extern int crack_mesh_line_width;
extern int crack_front_edge_line_width;
extern int segment_line_width;
extern int point_size;
extern int contour_type;
extern int reversed_contour_colors_flag;
extern int gray_contour_color_below_magenta_one_flag;
extern int num_contour_grades;
extern int contour_bar_flag;
extern char contour_bar_label[];
extern char contour_bar_label2[];
extern int auto_contour_bar_format_flag;
extern char contour_bar_format[];
extern double contour_bar_scale_factor;
extern int auto_contour_field_range_flag;
extern double contour_field_max;
extern double contour_field_min;
extern int integer_contour_field_mod;
extern int contour_field_id;
extern int displacement_field_id;
extern int point_mark;
extern int points_reversed_contour_colors_flag;
extern int points_contour_bar_flag;
extern char points_contour_bar_label[];
extern char points_contour_bar_label2[];
extern int points_auto_contour_bar_format_flag;
extern char points_contour_bar_format[];
extern double points_contour_bar_scale_factor;
extern int auto_points_field_range_flag;
extern double points_field_max;
extern double points_field_min;
extern int points_field_id;

/*
 * Set options
 */
void setOptions(int argc, char *argv[])
{
    /* Set options from arguments */
    setOptionsFromArguments(argc, argv);

    /* Input restart file */
    inputRestartFile();
}

/*
 * Set options from arguments
 */
void setOptionsFromArguments(int argc, char *argv[])
{
    char *endptr;
    int i;

    for (i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-h") == 0)
        {
            displayHelp();
            exit(EXIT_SUCCESS);
        }
        else if (strcmp(argv[i], "-e") == 0)
        {
            strcpy(filename_visible_elements, argv[++i]);
        }
        else if (strcmp(argv[i], "-E") == 0)
        {
            strcpy(filename_invisible_elements, argv[++i]);
        }
        else if (strcmp(argv[i], "-p") == 0)
        {
            strcpy(filename_points, argv[++i]);
        }
        else if (strcmp(argv[i], "-P") == 0)
        {
            strcpy(filename_points_fields, argv[++i]);
        }
        else if (strcmp(argv[i], "-r") == 0)
        {
            strcpy(filename_restart, argv[++i]);
        }
        else if (strcmp(argv[i], "-b") == 0)
        {
            batch_mode_flag = 1;
            strcpy(filename_view_image, argv[++i]);
        }
        else if (strcmp(argv[i], "-w") == 0)
        {
            window_size = (int)strtol(argv[++i], &endptr, 10);
            if (*endptr != '\0')
            {
                fprintf(stderr, "Error: invalid option -w %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-n") == 0)
        {
            division_number = (int)strtol(argv[++i], &endptr, 10);
            if (*endptr != '\0')
            {
                fprintf(stderr, "Error: invalid option -n %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-t") == 0)
        {
            first_time_step   = (int)strtol(argv[++i], &endptr, 10);
            if (*endptr != '\0')
            {
                fprintf(stderr, "Error: invalid option -t %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (strcmp(argv[i], "-T") == 0)
        {
            time_step_increment = (int)strtol(argv[++i], &endptr, 10);
            if (*endptr != '\0')
            {
                fprintf(stderr, "Error: invalid option -T %s\n", argv[i]);
                exit(EXIT_FAILURE);
            }
        }
        else if (argv[i][0] == '-')
        {
            fprintf(stderr, "Error: unknown option %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
        else
        {
            break;
        }
    }

    /* Display help if arguments not enough */
    if (i >= argc)
    {
        displayHelp();
        exit(EXIT_FAILURE);
    }

    /* Get filenames */
    filename_mesh = argv[i];
    num_filenames_fields = argc - (i + 1);
    filenames_fields = &argv[i + 1];
}

/*
 * Display help
 */
void displayHelp(void)
{
    printf("Usage:\n"
           "  zarusoba_viewer [options] mesh.adv [merged_advhddm_out.adv]\n"
           "  zarusoba_viewer [options] mesh.msh [hddmmrg.dat ...]\n"
           "  zarusoba_viewer [options] mesh.pch\n"
           "  zarusoba_viewer [options] mesh.pcm\n"
           "\n");

    printf("Input files:\n"
           "  mesh.adv     ADVENTURE-format mesh file\n"
           "  mesh.msh     Text-based mesh file\n"
           "  mesh.emsh    Text-based extended mesh file\n"
           "  mesh.pch     Text-based surface patch mesh file\n"
           "  mesh.pcm     Text-based multi-material surface patch mesh file\n"
           "  merged_advhddm_out.adv\n"
           "               ADVENTURE-format fields file generated by tools/merge_advhddm_out\n"
           "  hddmmrg.dat  Text-based fields file generated by hddmmrg in ADVENTURE_Solid\n"
           "\n");

    printf("Options:\n");
    printf("  -h                         Display help\n");
    printf("  -e visible_elements.txt    Input visible elements\n");
    printf("  -E invisible_elements.txt  Input invisible elements\n");
    printf("  -p points.txt              Input points\n");
    printf("  -P points_fields.txt       Input points fields\n");
    printf("  -r restart.txt             Input restart file\n");
    printf("  -b image.png               Enable batch mode and save view image\n");
    printf("  -w size                    Set window size [%d]\n", DEFAULT_WINDOW_SIZE);
    printf("  -n number                  Set division number to detect edges [%d]\n", DEFAULT_DIVISION_NUMBER);
    printf("  -t step                    Set first time step [%d]\n", DEFAULT_FIRST_TIME_STEP);
    printf("  -T step                    Set time step increment [%d]\n", DEFAULT_TIME_STEP_INCREMENT);
}

/*
 * Input restart file
 */
void inputRestartFile(void)
{
    FILE *fp;
    char buffer[FILE_INPUT_BUFFER_SIZE];
    char option_name[256], option_value[256];

    if (strlen(filename_restart) == 0)
    {
        strcpy(filename_restart, "restart.txt");
        return;
    }

    /* Zeroize flags */
    auto_view_position_flag = 0;
    auto_contour_field_range_flag = 0;
    auto_points_field_range_flag = 0;

    /* Open file */
    printf("Input %s\n", filename_restart);
    if ((fp = fopen(filename_restart, "r")) == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename_restart);
        exit(EXIT_FAILURE);
    }

    /* Read restart file */
    while (fgets(buffer, FILE_INPUT_BUFFER_SIZE, fp) != NULL)
    {
        switch (sscanf(buffer, "%s %[^\n]", option_name, option_value))
        {
        case 1:
            strcpy(option_value, "");
            break;
        case 2:
            break;
        default:
            fprintf(stderr, "Warning: failed to read restart options in %s\n", filename_restart);
        }

        if (strcmp(option_name, "view_size") == 0)
        {
            view_size = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "view_center_x") == 0)
        {
            view_center_x = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "view_center_y") == 0)
        {
            view_center_y = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "view_center_z") == 0)
        {
            view_center_z = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "view_direction_x") == 0)
        {
            view_direction_x = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "view_direction_y") == 0)
        {
            view_direction_y = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "view_direction_z") == 0)
        {
            view_direction_z = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "view_up_vector_x") == 0)
        {
            view_up_vector_x = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "view_up_vector_y") == 0)
        {
            view_up_vector_y = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "view_up_vector_z") == 0)
        {
            view_up_vector_z = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "section_a") == 0)
        {
            section_a = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "section_b") == 0)
        {
            section_b = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "section_c") == 0)
        {
            section_c = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "section_d") == 0)
        {
            section_d = atof(option_value);
            auto_view_position_flag++;
        }
        else if (strcmp(option_name, "plot_type") == 0)
            plot_type = atoi(option_value);
        else if (strcmp(option_name, "mesh_flag") == 0)
            mesh_flag = atoi(option_value);
        else if (strcmp(option_name, "edges_flag") == 0)
            edges_flag = atoi(option_value);
        else if (strcmp(option_name, "crack_mesh_flag") == 0)
            crack_mesh_flag = atoi(option_value);
        else if (strcmp(option_name, "crack_front_edges_flag") == 0)
            crack_front_edges_flag = atoi(option_value);
        else if (strcmp(option_name, "segments_flag") == 0)
            segments_flag = atoi(option_value);
        else if (strcmp(option_name, "deformed_mesh_flag") == 0)
            deformed_mesh_flag = atoi(option_value);
        else if (strcmp(option_name, "deformed_edges_flag") == 0)
            deformed_edges_flag = atoi(option_value);
        else if (strcmp(option_name, "deformed_crack_mesh_flag") == 0)
            deformed_crack_mesh_flag = atoi(option_value);
        else if (strcmp(option_name, "deformed_crack_front_edges_flag") == 0)
            deformed_crack_front_edges_flag = atoi(option_value);
        else if (strcmp(option_name, "deformed_segments_flag") == 0)
            deformed_segments_flag = atoi(option_value);
        else if (strcmp(option_name, "deformation_scale_factor") == 0)
            deformation_scale_factor = atof(option_value);
        else if (strcmp(option_name, "points_type") == 0)
            points_type = atoi(option_value);
        else if (strcmp(option_name, "title") == 0)
            strcpy(title, option_value);
        else if (strcmp(option_name, "coord_axes_flag") == 0)
            coord_axes_flag = atoi(option_value);
        else if (strcmp(option_name, "coord_axes_at_origin_flag") == 0)
            coord_axes_at_origin_flag = atoi(option_value);
        else if (strcmp(option_name, "other_mode") == 0)
            other_mode = atoi(option_value);
        else if (strcmp(option_name, "marked_node") == 0)
            marked_node = atoi(option_value);
        else if (strcmp(option_name, "marked_element") == 0)
            marked_element = atoi(option_value);
        else if (strcmp(option_name, "marked_point") == 0)
            marked_point = atoi(option_value);
        else if (strcmp(option_name, "scale_bar_flag") == 0)
            scale_bar_flag = atoi(option_value);
        else if (strcmp(option_name, "scale_bar_format") == 0)
            strcpy(scale_bar_format, option_value);
        else if (strcmp(option_name, "scale_bar_scale_factor") == 0)
            scale_bar_scale_factor = atof(option_value);
        else if (strcmp(option_name, "section_flag") == 0)
            section_flag = atoi(option_value);
        else if (strcmp(option_name, "section_type") == 0)
            section_type = atoi(option_value);
        else if (strcmp(option_name, "volume_color") == 0)
            strcpy(volume_color, option_value);
        else if (strcmp(option_name, "mesh_color") == 0)
            strcpy(mesh_color, option_value);
        else if (strcmp(option_name, "edge_color") == 0)
            strcpy(edge_color, option_value);
        else if (strcmp(option_name, "crack_mesh_color") == 0)
            strcpy(crack_mesh_color, option_value);
        else if (strcmp(option_name, "crack_front_edge_color") == 0)
            strcpy(crack_front_edge_color, option_value);
        else if (strcmp(option_name, "segment_color") == 0)
            strcpy(segment_color, option_value);
        else if (strcmp(option_name, "point_color") == 0)
            strcpy(point_color, option_value);
        else if (strcmp(option_name, "mesh_line_width") == 0)
            mesh_line_width = atoi(option_value);
        else if (strcmp(option_name, "edge_line_width") == 0)
            edge_line_width = atoi(option_value);
        else if (strcmp(option_name, "crack_mesh_line_width") == 0)
            crack_mesh_line_width = atoi(option_value);
        else if (strcmp(option_name, "crack_front_edge_line_width") == 0)
            crack_front_edge_line_width = atoi(option_value);
        else if (strcmp(option_name, "segment_line_width") == 0)
            segment_line_width = atoi(option_value);
        else if (strcmp(option_name, "point_size") == 0)
            point_size = atoi(option_value);
        else if (strcmp(option_name, "contour_type") == 0)
            contour_type = atoi(option_value);
        else if (strcmp(option_name, "reversed_contour_colors_flag") == 0)
            reversed_contour_colors_flag = atoi(option_value);
        else if (strcmp(option_name, "gray_contour_color_below_magenta_one_flag") == 0)
            gray_contour_color_below_magenta_one_flag = atoi(option_value);
        else if (strcmp(option_name, "num_contour_grades") == 0)
            num_contour_grades = atoi(option_value);
        else if (strcmp(option_name, "contour_bar_flag") == 0)
            contour_bar_flag = atoi(option_value);
        else if (strcmp(option_name, "contour_bar_label") == 0)
            strcpy(contour_bar_label2, option_value);
        else if (strcmp(option_name, "auto_contour_bar_format_flag") == 0)
            auto_contour_bar_format_flag = atoi(option_value);
        else if (strcmp(option_name, "contour_bar_format") == 0)
            strcpy(contour_bar_format, option_value);
        else if (strcmp(option_name, "contour_bar_scale_factor") == 0)
            contour_bar_scale_factor = atof(option_value);
        else if (strcmp(option_name, "contour_field_max") == 0)
        {
            contour_field_max = atof(option_value);
            auto_contour_field_range_flag++;
        }
        else if (strcmp(option_name, "contour_field_min") == 0)
        {
            contour_field_min = atof(option_value);
            auto_contour_field_range_flag++;
        }
        else if (strcmp(option_name, "integer_contour_field_mod") == 0)
            integer_contour_field_mod = atoi(option_value);
        else if (strcmp(option_name, "contour_field_id") == 0)
            contour_field_id = atoi(option_value);
        else if (strcmp(option_name, "displacement_field_id") == 0)
            displacement_field_id = atoi(option_value);
        else if (strcmp(option_name, "point_mark") == 0)
            point_mark = atoi(option_value);
        else if (strcmp(option_name, "points_reversed_contour_colors_flag") == 0)
            points_reversed_contour_colors_flag = atoi(option_value);
        else if (strcmp(option_name, "points_contour_bar_flag") == 0)
            points_contour_bar_flag = atoi(option_value);
        else if (strcmp(option_name, "points_contour_bar_label") == 0)
            strcpy(points_contour_bar_label2, option_value);
        else if (strcmp(option_name, "points_auto_contour_bar_format_flag") == 0)
            points_auto_contour_bar_format_flag = atoi(option_value);
        else if (strcmp(option_name, "points_contour_bar_format") == 0)
            strcpy(points_contour_bar_format, option_value);
        else if (strcmp(option_name, "points_contour_bar_scale_factor") == 0)
            points_contour_bar_scale_factor = atof(option_value);
        else if (strcmp(option_name, "points_field_max") == 0)
        {
            points_field_max = atof(option_value);
            auto_points_field_range_flag++;
        }
        else if (strcmp(option_name, "points_field_min") == 0)
        {
            points_field_min = atof(option_value);
            auto_points_field_range_flag++;
        }
        else if (strcmp(option_name, "points_field_id") == 0)
            points_field_id = atoi(option_value);
        else
            fprintf(stderr, "Warning: invalid restart option %s\n", option_name);
    }

    /* Set flags */
    auto_view_position_flag = (auto_view_position_flag == 14)? 0: 1;
    auto_contour_field_range_flag = (auto_contour_field_range_flag == 2)? 0: 1;
    auto_points_field_range_flag = (auto_points_field_range_flag == 2)? 0: 1;

    /* Close file */
    fclose(fp);
}

/*
 * Output restart file
 */
void outputRestartFile(void)
{
    FILE *fp;

    /* Open file */
    if ((fp = fopen(filename_restart, "w")) == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", filename_restart);
        exit(EXIT_FAILURE);
    }

    /* Write restart options */
    fprintf(fp, "%s %.15E\n", "view_size", view_size);
    fprintf(fp, "%s %.15E\n", "view_center_x", view_center_x);
    fprintf(fp, "%s %.15E\n", "view_center_y", view_center_y);
    fprintf(fp, "%s %.15E\n", "view_center_z", view_center_z);
    fprintf(fp, "%s %.15E\n", "view_direction_x", view_direction_x);
    fprintf(fp, "%s %.15E\n", "view_direction_y", view_direction_y);
    fprintf(fp, "%s %.15E\n", "view_direction_z", view_direction_z);
    fprintf(fp, "%s %.15E\n", "view_up_vector_x", view_up_vector_x);
    fprintf(fp, "%s %.15E\n", "view_up_vector_y", view_up_vector_y);
    fprintf(fp, "%s %.15E\n", "view_up_vector_z", view_up_vector_z);
    fprintf(fp, "%s %.15E\n", "section_a", section_a);
    fprintf(fp, "%s %.15E\n", "section_b", section_b);
    fprintf(fp, "%s %.15E\n", "section_c", section_c);
    fprintf(fp, "%s %.15E\n", "section_d", section_d);
    fprintf(fp, "%s %d\n", "plot_type", plot_type);
    fprintf(fp, "%s %d\n", "mesh_flag", mesh_flag);
    fprintf(fp, "%s %d\n", "edges_flag", edges_flag);
    fprintf(fp, "%s %d\n", "crack_mesh_flag", crack_mesh_flag);
    fprintf(fp, "%s %d\n", "crack_front_edges_flag", crack_front_edges_flag);
    fprintf(fp, "%s %d\n", "segments_flag", segments_flag);
    fprintf(fp, "%s %d\n", "deformed_mesh_flag", deformed_mesh_flag);
    fprintf(fp, "%s %d\n", "deformed_edges_flag", deformed_edges_flag);
    fprintf(fp, "%s %d\n", "deformed_crack_mesh_flag", deformed_crack_mesh_flag);
    fprintf(fp, "%s %d\n", "deformed_crack_front_edges_flag", deformed_crack_front_edges_flag);
    fprintf(fp, "%s %d\n", "deformed_segments_flag", deformed_segments_flag);
    fprintf(fp, "%s %.15E\n", "deformation_scale_factor", deformation_scale_factor);
    fprintf(fp, "%s %d\n", "points_type", points_type);
    fprintf(fp, "%s %s\n", "title", title);
    fprintf(fp, "%s %d\n", "coord_axes_flag", coord_axes_flag);
    fprintf(fp, "%s %d\n", "coord_axes_at_origin_flag", coord_axes_at_origin_flag);
    fprintf(fp, "%s %d\n", "other_mode", other_mode);
    fprintf(fp, "%s %d\n", "marked_node", marked_node);
    fprintf(fp, "%s %d\n", "marked_element", marked_element);
    fprintf(fp, "%s %d\n", "marked_point", marked_point);
    fprintf(fp, "%s %d\n", "scale_bar_flag", scale_bar_flag);
    fprintf(fp, "%s %s\n", "scale_bar_format", scale_bar_format);
    fprintf(fp, "%s %.15E\n", "scale_bar_scale_factor", scale_bar_scale_factor);
    fprintf(fp, "%s %d\n", "section_flag", section_flag);
    fprintf(fp, "%s %d\n", "section_type", section_type);
    fprintf(fp, "%s %s\n", "volume_color", volume_color);
    fprintf(fp, "%s %s\n", "mesh_color", mesh_color);
    fprintf(fp, "%s %s\n", "edge_color", edge_color);
    fprintf(fp, "%s %s\n", "crack_mesh_color", crack_mesh_color);
    fprintf(fp, "%s %s\n", "crack_front_edge_color", crack_front_edge_color);
    fprintf(fp, "%s %s\n", "segment_color", segment_color);
    fprintf(fp, "%s %s\n", "point_color", point_color);
    fprintf(fp, "%s %d\n", "mesh_line_width", mesh_line_width);
    fprintf(fp, "%s %d\n", "edge_line_width", edge_line_width);
    fprintf(fp, "%s %d\n", "crack_mesh_line_width", crack_mesh_line_width);
    fprintf(fp, "%s %d\n", "crack_front_edge_line_width", crack_front_edge_line_width);
    fprintf(fp, "%s %d\n", "segment_line_width", segment_line_width);
    fprintf(fp, "%s %d\n", "point_size", point_size);
    fprintf(fp, "%s %d\n", "contour_type", contour_type);
    fprintf(fp, "%s %d\n", "reversed_contour_colors_flag", reversed_contour_colors_flag);
    fprintf(fp, "%s %d\n", "gray_contour_color_below_magenta_one_flag", gray_contour_color_below_magenta_one_flag);
    fprintf(fp, "%s %d\n", "num_contour_grades", num_contour_grades);
    fprintf(fp, "%s %d\n", "contour_bar_flag", contour_bar_flag);
    fprintf(fp, "%s %s\n", "contour_bar_label", contour_bar_label);
    fprintf(fp, "%s %d\n", "auto_contour_bar_format_flag", auto_contour_bar_format_flag);
    fprintf(fp, "%s %s\n", "contour_bar_format", contour_bar_format);
    fprintf(fp, "%s %.15E\n", "contour_bar_scale_factor", contour_bar_scale_factor);
    fprintf(fp, "%s %.15E\n", "contour_field_max", contour_field_max);
    fprintf(fp, "%s %.15E\n", "contour_field_min", contour_field_min);
    fprintf(fp, "%s %d\n", "integer_contour_field_mod", integer_contour_field_mod);
    fprintf(fp, "%s %d\n", "contour_field_id", contour_field_id);
    fprintf(fp, "%s %d\n", "displacement_field_id", displacement_field_id);
    fprintf(fp, "%s %d\n", "point_mark", point_mark);
    fprintf(fp, "%s %d\n", "points_reversed_contour_colors_flag", points_reversed_contour_colors_flag);
    fprintf(fp, "%s %d\n", "points_contour_bar_flag", points_contour_bar_flag);
    fprintf(fp, "%s %s\n", "points_contour_bar_label", points_contour_bar_label);
    fprintf(fp, "%s %d\n", "points_auto_contour_bar_format_flag", points_auto_contour_bar_format_flag);
    fprintf(fp, "%s %s\n", "points_contour_bar_format", points_contour_bar_format);
    fprintf(fp, "%s %.15E\n", "points_contour_bar_scale_factor", points_contour_bar_scale_factor);
    fprintf(fp, "%s %.15E\n", "points_field_max", points_field_max);
    fprintf(fp, "%s %.15E\n", "points_field_min", points_field_min);
    fprintf(fp, "%s %d\n", "points_field_id", points_field_id);

    /* Close file */
    fclose(fp);
}
