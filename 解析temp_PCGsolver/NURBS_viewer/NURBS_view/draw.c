#include <stdio.h>
#include <string.h>
#include <math.h>
#include <autogl.h>

#include "draw.h"

#include "constant.h"

extern double section_a;
extern double section_b;
extern double section_c;
extern double section_d;

extern int section_type;
extern int point_size;
extern int contour_type;
extern int num_contour_grades;
extern double contour_field_max;
extern double contour_field_min;
extern int point_mark;

static void drawExactContourTetrahedronSection(
	double x0, double y0, double z0, double value0,
	double x1, double y1, double z1, double value1,
	double x2, double y2, double z2, double value2,
	double x3, double y3, double z3, double value3);
static void drawElementalContourTetrahedronSection(
	double x0, double y0, double z0, double value0,
	double x1, double y1, double z1, double value1,
	double x2, double y2, double z2, double value2,
	double x3, double y3, double z3, double value3);
static void drawExactContourHexahedronSection(
	double x0, double y0, double z0, double value0,
	double x1, double y1, double z1, double value1,
	double x2, double y2, double z2, double value2,
	double x3, double y3, double z3, double value3,
	double x4, double y4, double z4, double value4,
	double x5, double y5, double z5, double value5,
	double x6, double y6, double z6, double value6,
	double x7, double y7, double z7, double value7);
static void drawElementalContourHexahedronSection(
	double x0, double y0, double z0, double value0,
	double x1, double y1, double z1, double value1,
	double x2, double y2, double z2, double value2,
	double x3, double y3, double z3, double value3,
	double x4, double y4, double z4, double value4,
	double x5, double y5, double z5, double value5,
	double x6, double y6, double z6, double value6,
	double x7, double y7, double z7, double value7);

/*
 * Draw line
 */
void drawLine(double x0, double y0, double z0,
			  double x1, double y1, double z1) {
	AutoGL_DrawLine(COORDS_SCALE_FACTOR * x0,
					COORDS_SCALE_FACTOR * y0,
					COORDS_SCALE_FACTOR * z0,
					COORDS_SCALE_FACTOR * x1,
					COORDS_SCALE_FACTOR * y1,
					COORDS_SCALE_FACTOR * z1);
}

/*
 * Draw sectioned line
 */
void drawSectionedLine(double x0, double y0, double z0,
					   double x1, double y1, double z1) {
	double h0 = section_a * x0 + section_b * y0 + section_c * z0 + section_d;
	double h1 = section_a * x1 + section_b * y1 + section_c * z1 + section_d;

	int flag0 = 0;
	int flag1 = 0;

	double xx0 = x0;
	double yy0 = y0;
	double zz0 = z0;
	double hh0 = h0;
	double xx1 = x1;
	double yy1 = y1;
	double zz1 = z1;
	double hh1 = h1;

	switch (section_type) {
	case LESS_EQUAL_SECTION_TYPE:
		flag0 = (h0 < 0.0);
		flag1 = (h1 < 0.0);
		break;
	case GREATER_EQUAL_SECTION_TYPE:
		flag0 = (h0 > 0.0);
		flag1 = (h1 > 0.0);
		break;
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
		flag0 = (h0 <= 0.0);
		flag1 = (h1 <= 0.0);
		break;
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		flag0 = (h0 >= 0.0);
		flag1 = (h1 >= 0.0);
		break;
	default:
		return;
	}

	if (!flag0 && !flag1) {
		return;
	} else if (flag0 && flag1) {
		drawLine(x0, y0, z0,
				 x1, y1, z1);
		return;
	}

	switch (section_type) {
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		if (!flag0 || !flag1)
			return;
	}

	h0 = fabs(h0);
	h1 = fabs(h1);

	hh0 = h0;
	hh1 = h1;

	if (!flag0 && flag1) {
		flag0 = 1; flag1 = 0;
		xx0 = x1; yy0 = y1; zz0 = z1; hh0 = h1;
		xx1 = x0; yy1 = y0; zz1 = z0; hh1 = h0;
	}

	drawLine(xx0,
			 yy0,
			 zz0,
			 (hh0 * xx1 + hh1 * xx0) / (hh0 + hh1),
			 (hh0 * yy1 + hh1 * yy0) / (hh0 + hh1),
			 (hh0 * zz1 + hh1 * zz0) / (hh0 + hh1));
}

/*
 * Draw tetrahedron section lines
 */
void drawTetrahedronSectionLines(double x0, double y0, double z0,
								 double x1, double y1, double z1,
								 double x2, double y2, double z2,
								 double x3, double y3, double z3) {
	double h0 = section_a * x0 + section_b * y0 + section_c * z0 + section_d;
	double h1 = section_a * x1 + section_b * y1 + section_c * z1 + section_d;
	double h2 = section_a * x2 + section_b * y2 + section_c * z2 + section_d;
	double h3 = section_a * x3 + section_b * y3 + section_c * z3 + section_d;

	int flag0 = 0;
	int flag1 = 0;
	int flag2 = 0;
	int flag3 = 0;

	switch (section_type) {
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
		flag0 = (h0 <= 0.0);
		flag1 = (h1 <= 0.0);
		flag2 = (h2 <= 0.0);
		flag3 = (h3 <= 0.0);
		break;
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		flag0 = (h0 >= 0.0);
		flag1 = (h1 >= 0.0);
		flag2 = (h2 >= 0.0);
		flag3 = (h3 >= 0.0);
		break;
	default:
		return;
	}

	if (flag0 && flag1 && flag2 && flag3)
		return;

	if (!flag0 && !flag1 && !flag2 && !flag3)
		return;

	if (flag1 && flag2 && flag3) {
		drawLine(x1, y1, z1,
				 x2, y2, z2);
		drawLine(x2, y2, z2,
				 x3, y3, z3);
		drawLine(x3, y3, z3,
				 x1, y1, z1);
	}

	if (flag0 && flag3 && flag2) {
		drawLine(x0, y0, z0,
				 x3, y3, z3);
		drawLine(x3, y3, z3,
				 x2, y2, z2);
		drawLine(x2, y2, z2,
				 x0, y0, z0);
	}

	if (flag0 && flag1 && flag3) {
		drawLine(x0, y0, z0,
				 x1, y1, z1);
		drawLine(x1, y1, z1,
				 x3, y3, z3);
		drawLine(x3, y3, z3,
				 x0, y0, z0);
	}

	if (flag0 && flag2 && flag1) {
		drawLine(x0, y0, z0,
				 x2, y2, z2);
		drawLine(x2, y2, z2,
				 x1, y1, z1);
		drawLine(x1, y1, z1,
				 x0, y0, z0);
	}
}

/*
 * Draw hexahedron section lines
 */
void drawHexahedronSectionLines(double x0, double y0, double z0,
								double x1, double y1, double z1,
								double x2, double y2, double z2,
								double x3, double y3, double z3,
								double x4, double y4, double z4,
								double x5, double y5, double z5,
								double x6, double y6, double z6,
								double x7, double y7, double z7) {
	double h0 = section_a * x0 + section_b * y0 + section_c * z0 + section_d;
	double h1 = section_a * x1 + section_b * y1 + section_c * z1 + section_d;
	double h2 = section_a * x2 + section_b * y2 + section_c * z2 + section_d;
	double h3 = section_a * x3 + section_b * y3 + section_c * z3 + section_d;
	double h4 = section_a * x4 + section_b * y4 + section_c * z4 + section_d;
	double h5 = section_a * x5 + section_b * y5 + section_c * z5 + section_d;
	double h6 = section_a * x6 + section_b * y6 + section_c * z6 + section_d;
	double h7 = section_a * x7 + section_b * y7 + section_c * z7 + section_d;

	int flag0 = 0;
	int flag1 = 0;
	int flag2 = 0;
	int flag3 = 0;
	int flag4 = 0;
	int flag5 = 0;
	int flag6 = 0;
	int flag7 = 0;

	switch (section_type) {
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
		flag0 = (h0 <= 0.0);
		flag1 = (h1 <= 0.0);
		flag2 = (h2 <= 0.0);
		flag3 = (h3 <= 0.0);
		flag4 = (h4 <= 0.0);
		flag5 = (h5 <= 0.0);
		flag6 = (h6 <= 0.0);
		flag7 = (h7 <= 0.0);
		break;
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		flag0 = (h0 >= 0.0);
		flag1 = (h1 >= 0.0);
		flag2 = (h2 >= 0.0);
		flag3 = (h3 >= 0.0);
		flag4 = (h4 >= 0.0);
		flag5 = (h5 >= 0.0);
		flag6 = (h6 >= 0.0);
		flag7 = (h7 >= 0.0);
		break;
	default:
		return;
	}

	if (flag0 && flag1 && flag2 && flag3 && flag4 && flag5 && flag6 && flag7)
		return;

	if (!flag0 && !flag1 && !flag2 && !flag3 && !flag4 && !flag5 && !flag6 && !flag7)
		return;

	if (flag0 && flag3 && flag2 && flag1) {
		drawLine(x0, y0, z0,
				 x3, y3, z3);
		drawLine(x3, y3, z3,
				 x2, y2, z2);
		drawLine(x2, y2, z2,
				 x1, y1, z1);
		drawLine(x1, y1, z1,
				 x0, y0, z0);
	}

	if (flag4 && flag5 && flag6 && flag7) {
		drawLine(x4, y4, z4,
				 x5, y5, z5);
		drawLine(x5, y5, z5,
				 x6, y6, z6);
		drawLine(x6, y6, z6,
				 x7, y7, z7);
		drawLine(x7, y7, z7,
				 x4, y4, z4);
	}

	if (flag0 && flag1 && flag5 && flag4) {
		drawLine(x0, y0, z0,
				 x1, y1, z1);
		drawLine(x1, y1, z1,
				 x5, y5, z5);
		drawLine(x5, y5, z5,
				 x4, y4, z4);
		drawLine(x4, y4, z4,
				 x0, y0, z0);
	}

	if (flag1 && flag2 && flag6 && flag5) {
		drawLine(x1, y1, z1,
				 x2, y2, z2);
		drawLine(x2, y2, z2,
				 x6, y6, z6);
		drawLine(x6, y6, z6,
				 x5, y5, z5);
		drawLine(x5, y5, z5,
				 x1, y1, z1);
	}

	if (flag2 && flag3 && flag7 && flag6) {
		drawLine(x2, y2, z2,
				 x3, y3, z3);
		drawLine(x3, y3, z3,
				 x7, y7, z7);
		drawLine(x7, y7, z7,
				 x6, y6, z6);
		drawLine(x6, y6, z6,
				 x2, y2, z2);
	}

	if (flag0 && flag4 && flag7 && flag3) {
		drawLine(x0, y0, z0,
				 x4, y4, z4);
		drawLine(x4, y4, z4,
				 x7, y7, z7);
		drawLine(x7, y7, z7,
				 x3, y3, z3);
		drawLine(x3, y3, z3,
				 x0, y0, z0);
	}
}

/*
 * Draw triangle
 */
void drawTriangle(double x0, double y0, double z0,
				  double x1, double y1, double z1,
				  double x2, double y2, double z2) {
	AutoGL_DrawTriangle(COORDS_SCALE_FACTOR * x0,
						COORDS_SCALE_FACTOR * y0,
						COORDS_SCALE_FACTOR * z0,
						COORDS_SCALE_FACTOR * x1,
						COORDS_SCALE_FACTOR * y1,
						COORDS_SCALE_FACTOR * z1,
						COORDS_SCALE_FACTOR * x2,
						COORDS_SCALE_FACTOR * y2,
						COORDS_SCALE_FACTOR * z2);
}

/*
 * Draw quadrangle
 */
void drawQuadrangle(double x0, double y0, double z0,
					double x1, double y1, double z1,
					double x2, double y2, double z2,
					double x3, double y3, double z3) {
	AutoGL_DrawQuadrangle(COORDS_SCALE_FACTOR * x0,
						  COORDS_SCALE_FACTOR * y0,
						  COORDS_SCALE_FACTOR * z0,
						  COORDS_SCALE_FACTOR * x1,
						  COORDS_SCALE_FACTOR * y1,
						  COORDS_SCALE_FACTOR * z1,
						  COORDS_SCALE_FACTOR * x2,
						  COORDS_SCALE_FACTOR * y2,
						  COORDS_SCALE_FACTOR * z2,
						  COORDS_SCALE_FACTOR * x3,
						  COORDS_SCALE_FACTOR * y3,
						  COORDS_SCALE_FACTOR * z3);
}

/*
 * Draw sectioned triangle
 */
void drawSectionedTriangle(double x0, double y0, double z0,
						   double x1, double y1, double z1,
						   double x2, double y2, double z2) {
	double h0 = section_a * x0 + section_b * y0 + section_c * z0 + section_d;
	double h1 = section_a * x1 + section_b * y1 + section_c * z1 + section_d;
	double h2 = section_a * x2 + section_b * y2 + section_c * z2 + section_d;

	int flag0 = 0;
	int flag1 = 0;
	int flag2 = 0;

	double xx0 = x0;
	double yy0 = y0;
	double zz0 = z0;
	double hh0 = h0;
	double xx1 = x1;
	double yy1 = y1;
	double zz1 = z1;
	double hh1 = h1;
	double xx2 = x2;
	double yy2 = y2;
	double zz2 = z2;
	double hh2 = h2;

	switch (section_type) {
	case LESS_EQUAL_SECTION_TYPE:
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
		flag0 = (h0 <= 0.0);
		flag1 = (h1 <= 0.0);
		flag2 = (h2 <= 0.0);
		break;
	case GREATER_EQUAL_SECTION_TYPE:
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		flag0 = (h0 >= 0.0);
		flag1 = (h1 >= 0.0);
		flag2 = (h2 >= 0.0);
		break;
	default:
		return;
	}

	if (!flag0 && !flag1 && !flag2) {
		return;
	} else if (flag0 && flag1 && flag2) {
		drawTriangle(x0, y0, z0,
					 x1, y1, z1,
					 x2, y2, z2);
		return;
	}

	switch (section_type) {
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		if (!flag0 || !flag1 || !flag2)
			return;
	}

	h0 = fabs(h0);
	h1 = fabs(h1);
	h2 = fabs(h2);

	hh0 = h0;
	hh1 = h1;
	hh2 = h2;

	if ((!flag0 && flag1 && !flag2) || (!flag0 && flag1 && flag2)) {
		flag1 = flag2;
		flag0 = 1; flag2 = 0;
		xx0 = x1; yy0 = y1; zz0 = z1; hh0 = h1;
		xx1 = x2; yy1 = y2; zz1 = z2; hh1 = h2;
		xx2 = x0; yy2 = y0; zz2 = z0; hh2 = h0;

	} else if ((!flag0 && !flag1 && flag2) || (flag0 && !flag1 && flag2)) {
		flag1 = flag0;
		flag0 = 1; flag2 = 0;
		xx0 = x2; yy0 = y2; zz0 = z2; hh0 = h2;
		xx1 = x0; yy1 = y0; zz1 = z0; hh1 = h0;
		xx2 = x1; yy2 = y1; zz2 = z1; hh2 = h1;
	}

	if (flag0 && !flag1 && !flag2)
		drawTriangle(xx0,
					 yy0,
					 zz0,
					 (hh0 * xx1 + hh1 * xx0) / (hh0 + hh1),
					 (hh0 * yy1 + hh1 * yy0) / (hh0 + hh1),
					 (hh0 * zz1 + hh1 * zz0) / (hh0 + hh1),
					 (hh0 * xx2 + hh2 * xx0) / (hh0 + hh2),
					 (hh0 * yy2 + hh2 * yy0) / (hh0 + hh2),
					 (hh0 * zz2 + hh2 * zz0) / (hh0 + hh2));
	else if (flag0 && flag1 && !flag2)
		drawQuadrangle(xx0,
					   yy0,
					   zz0,
					   xx1,
					   yy1,
					   zz1,
					   (hh1 * xx2 + hh2 * xx1) / (hh1 + hh2),
					   (hh1 * yy2 + hh2 * yy1) / (hh1 + hh2),
					   (hh1 * zz2 + hh2 * zz1) / (hh1 + hh2),
					   (hh0 * xx2 + hh2 * xx0) / (hh0 + hh2),
					   (hh0 * yy2 + hh2 * yy0) / (hh0 + hh2),
					   (hh0 * zz2 + hh2 * zz0) / (hh0 + hh2));
}

/*
 * Draw sectioned quadrangle
 */
void drawSectionedQuadrangle(double x0, double y0, double z0,
							 double x1, double y1, double z1,
							 double x2, double y2, double z2,
							 double x3, double y3, double z3) {
	double h0 = section_a * x0 + section_b * y0 + section_c * z0 + section_d;
	double h1 = section_a * x1 + section_b * y1 + section_c * z1 + section_d;
	double h2 = section_a * x2 + section_b * y2 + section_c * z2 + section_d;
	double h3 = section_a * x3 + section_b * y3 + section_c * z3 + section_d;

	int flag0 = 0;
	int flag1 = 0;
	int flag2 = 0;
	int flag3 = 0;

	switch (section_type) {
	case LESS_EQUAL_SECTION_TYPE:
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
		flag0 = (h0 <= 0.0);
		flag1 = (h1 <= 0.0);
		flag2 = (h2 <= 0.0);
		flag3 = (h3 <= 0.0);
		break;
	case GREATER_EQUAL_SECTION_TYPE:
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		flag0 = (h0 >= 0.0);
		flag1 = (h1 >= 0.0);
		flag2 = (h2 >= 0.0);
		flag3 = (h3 >= 0.0);
		break;
	default:
		return;
	}

	if (!flag0 && !flag1 && !flag2 && !flag3) {
		return;
	} else if (flag0 && flag1 && flag2 && flag3) {
		drawQuadrangle(x0, y0, z0,
					   x1, y1, z1,
					   x2, y2, z2,
					   x3, y3, z3);
		return;
	}

	switch (section_type) {
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		if (!flag0 || !flag1 || !flag2 || !flag3)
			return;
	}

	drawSectionedTriangle(x0, y0, z0,
						  x1, y1, z1,
						  x2, y2, z2);
	drawSectionedTriangle(x2, y2, z2,
						  x3, y3, z3,
						  x0, y0, z0);
}

/*
 * Draw contour triangle
 */
void drawContourTriangle(double x0, double y0, double z0, double value0,
						 double x1, double y1, double z1, double value1,
						 double x2, double y2, double z2, double value2) {
	if (contour_type == SMOOTH_CONTOUR_TYPE)
		AutoGL_DrawSmoothContourTriangle(contour_field_min, contour_field_max,
										 COORDS_SCALE_FACTOR * x0,
										 COORDS_SCALE_FACTOR * y0,
										 COORDS_SCALE_FACTOR * z0,
										 value0,
										 COORDS_SCALE_FACTOR * x1,
										 COORDS_SCALE_FACTOR * y1,
										 COORDS_SCALE_FACTOR * z1,
										 value1,
										 COORDS_SCALE_FACTOR * x2,
										 COORDS_SCALE_FACTOR * y2,
										 COORDS_SCALE_FACTOR * z2,
										 value2);
	else if (contour_type == BAND_CONTOUR_TYPE)
		AutoGL_DrawBandContourTriangle(contour_field_min, contour_field_max,
									   num_contour_grades,
									   COORDS_SCALE_FACTOR * x0,
									   COORDS_SCALE_FACTOR * y0,
									   COORDS_SCALE_FACTOR * z0,
									   value0,
									   COORDS_SCALE_FACTOR * x1,
									   COORDS_SCALE_FACTOR * y1,
									   COORDS_SCALE_FACTOR * z1,
									   value1,
									   COORDS_SCALE_FACTOR * x2,
									   COORDS_SCALE_FACTOR * y2,
									   COORDS_SCALE_FACTOR * z2,
									   value2);
}

/*
 * Draw contour quadrangle
 */
void drawContourQuadrangle(double x0, double y0, double z0, double value0,
						   double x1, double y1, double z1, double value1,
						   double x2, double y2, double z2, double value2,
						   double x3, double y3, double z3, double value3) {
	if (contour_type == SMOOTH_CONTOUR_TYPE)
		AutoGL_DrawSmoothContourQuadrangle(contour_field_min, contour_field_max,
										   COORDS_SCALE_FACTOR * x0,
										   COORDS_SCALE_FACTOR * y0,
										   COORDS_SCALE_FACTOR * z0,
										   value0,
										   COORDS_SCALE_FACTOR * x1,
										   COORDS_SCALE_FACTOR * y1,
										   COORDS_SCALE_FACTOR * z1,
										   value1,
										   COORDS_SCALE_FACTOR * x2,
										   COORDS_SCALE_FACTOR * y2,
										   COORDS_SCALE_FACTOR * z2,
										   value2,
										   COORDS_SCALE_FACTOR * x3,
										   COORDS_SCALE_FACTOR * y3,
										   COORDS_SCALE_FACTOR * z3,
										   value3);
	else if (contour_type == BAND_CONTOUR_TYPE)
		AutoGL_DrawBandContourQuadrangle(contour_field_min, contour_field_max,
										 num_contour_grades,
										 COORDS_SCALE_FACTOR * x0,
										 COORDS_SCALE_FACTOR * y0,
										 COORDS_SCALE_FACTOR * z0,
										 value0,
										 COORDS_SCALE_FACTOR * x1,
										 COORDS_SCALE_FACTOR * y1,
										 COORDS_SCALE_FACTOR * z1,
										 value1,
										 COORDS_SCALE_FACTOR * x2,
										 COORDS_SCALE_FACTOR * y2,
										 COORDS_SCALE_FACTOR * z2,
										 value2,
										 COORDS_SCALE_FACTOR * x3,
										 COORDS_SCALE_FACTOR * y3,
										 COORDS_SCALE_FACTOR * z3,
										 value3);
}

/*
 * Draw sectioned contour triangle
 */
void drawSectionedContourTriangle(double x0, double y0, double z0, double value0,
								  double x1, double y1, double z1, double value1,
								  double x2, double y2, double z2, double value2) {
	double h0 = section_a * x0 + section_b * y0 + section_c * z0 + section_d;
	double h1 = section_a * x1 + section_b * y1 + section_c * z1 + section_d;
	double h2 = section_a * x2 + section_b * y2 + section_c * z2 + section_d;

	int flag0 = 0;
	int flag1 = 0;
	int flag2 = 0;

	double xx0 = x0;
	double yy0 = y0;
	double zz0 = z0;
	double vv0 = value0;
	double hh0 = h0;
	double xx1 = x1;
	double yy1 = y1;
	double zz1 = z1;
	double vv1 = value1;
	double hh1 = h1;
	double xx2 = x2;
	double yy2 = y2;
	double zz2 = z2;
	double vv2 = value2;
	double hh2 = h2;

	switch (section_type) {
	case LESS_EQUAL_SECTION_TYPE:
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
		flag0 = (h0 <= 0.0);
		flag1 = (h1 <= 0.0);
		flag2 = (h2 <= 0.0);
		break;
	case GREATER_EQUAL_SECTION_TYPE:
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		flag0 = (h0 >= 0.0);
		flag1 = (h1 >= 0.0);
		flag2 = (h2 >= 0.0);
		break;
	default:
		return;
	}

	if (!flag0 && !flag1 && !flag2) {
		return;
	} else if (flag0 && flag1 && flag2) {
		drawContourTriangle(x0, y0, z0, value0,
							x1, y1, z1, value1,
							x2, y2, z2, value2);
		return;
	}

	switch (section_type) {
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		if (!flag0 || !flag1 || !flag2)
			return;
	}

	h0 = fabs(h0);
	h1 = fabs(h1);
	h2 = fabs(h2);

	hh0 = h0;
	hh1 = h1;
	hh2 = h2;

	if ((!flag0 && flag1 && !flag2) || (!flag0 && flag1 && flag2)) {
		flag1 = flag2;
		flag0 = 1; flag2 = 0;
		xx0 = x1; yy0 = y1; zz0 = z1; vv0 = value1; hh0 = h1;
		xx1 = x2; yy1 = y2; zz1 = z2; vv1 = value2; hh1 = h2;
		xx2 = x0; yy2 = y0; zz2 = z0; vv2 = value0; hh2 = h0;

	} else if ((!flag0 && !flag1 && flag2) || (flag0 && !flag1 && flag2)) {
		flag1 = flag0;
		flag0 = 1; flag2 = 0;
		xx0 = x2; yy0 = y2; zz0 = z2; vv0 = value2; hh0 = h2;
		xx1 = x0; yy1 = y0; zz1 = z0; vv1 = value0; hh1 = h0;
		xx2 = x1; yy2 = y1; zz2 = z1; vv2 = value1; hh2 = h1;
	}

	if (flag0 && !flag1 && !flag2)
		drawContourTriangle(xx0,
							yy0,
							zz0,
							vv0,
							(hh0 * xx1 + hh1 * xx0) / (hh0 + hh1),
							(hh0 * yy1 + hh1 * yy0) / (hh0 + hh1),
							(hh0 * zz1 + hh1 * zz0) / (hh0 + hh1),
							(hh0 * vv1 + hh1 * vv0) / (hh0 + hh1),
							(hh0 * xx2 + hh2 * xx0) / (hh0 + hh2),
							(hh0 * yy2 + hh2 * yy0) / (hh0 + hh2),
							(hh0 * zz2 + hh2 * zz0) / (hh0 + hh2),
							(hh0 * vv2 + hh2 * vv0) / (hh0 + hh2));
	else if (flag0 && flag1 && !flag2)
		drawContourQuadrangle(xx0,
							  yy0,
							  zz0,
							  vv0,
							  xx1,
							  yy1,
							  zz1,
							  vv1,
							  (hh1 * xx2 + hh2 * xx1) / (hh1 + hh2),
							  (hh1 * yy2 + hh2 * yy1) / (hh1 + hh2),
							  (hh1 * zz2 + hh2 * zz1) / (hh1 + hh2),
							  (hh1 * vv2 + hh2 * vv1) / (hh1 + hh2),
							  (hh0 * xx2 + hh2 * xx0) / (hh0 + hh2),
							  (hh0 * yy2 + hh2 * yy0) / (hh0 + hh2),
							  (hh0 * zz2 + hh2 * zz0) / (hh0 + hh2),
							  (hh0 * vv2 + hh2 * vv0) / (hh0 + hh2));
}

/*
 * Draw sectioned contour quadrangle
 */
void drawSectionedContourQuadrangle(double x0, double y0, double z0, double value0,
									double x1, double y1, double z1, double value1,
									double x2, double y2, double z2, double value2,
									double x3, double y3, double z3, double value3) {
	double h0 = section_a * x0 + section_b * y0 + section_c * z0 + section_d;
	double h1 = section_a * x1 + section_b * y1 + section_c * z1 + section_d;
	double h2 = section_a * x2 + section_b * y2 + section_c * z2 + section_d;
	double h3 = section_a * x3 + section_b * y3 + section_c * z3 + section_d;

	int flag0 = 0;
	int flag1 = 0;
	int flag2 = 0;
	int flag3 = 0;

	switch (section_type) {
	case LESS_EQUAL_SECTION_TYPE:
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
		flag0 = (h0 <= 0.0);
		flag1 = (h1 <= 0.0);
		flag2 = (h2 <= 0.0);
		flag3 = (h3 <= 0.0);
		break;
	case GREATER_EQUAL_SECTION_TYPE:
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		flag0 = (h0 >= 0.0);
		flag1 = (h1 >= 0.0);
		flag2 = (h2 >= 0.0);
		flag3 = (h3 >= 0.0);
		break;
	default:
		return;
	}

	if (!flag0 && !flag1 && !flag2 && !flag3) {
		return;
	} else if (flag0 && flag1 && flag2 && flag3) {
		drawContourQuadrangle(x0, y0, z0, value0,
							  x1, y1, z1, value1,
							  x2, y2, z2, value2,
							  x3, y3, z3, value3);
		return;
	}

	switch (section_type) {
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		if (!flag0 || !flag1 || !flag2 || !flag3)
			return;
	}

	drawSectionedContourTriangle(x0, y0, z0, value0,
								 x1, y1, z1, value1,
								 x2, y2, z2, value2);
	drawSectionedContourTriangle(x2, y2, z2, value2,
								 x3, y3, z3, value3,
								 x0, y0, z0, value0);
}

/*
 * Draw contour tetrahedron section
 */
void drawContourTetrahedronSection(double x0, double y0, double z0, double value0,
								   double x1, double y1, double z1, double value1,
								   double x2, double y2, double z2, double value2,
								   double x3, double y3, double z3, double value3) {
	switch (section_type) {
	case EQUAL_SECTION_TYPE:
	case LESS_EQUAL_SECTION_TYPE:
	case GREATER_EQUAL_SECTION_TYPE:
		drawExactContourTetrahedronSection(x0, y0, z0, value0,
										   x1, y1, z1, value1,
										   x2, y2, z2, value2,
										   x3, y3, z3, value3);
		break;
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		drawElementalContourTetrahedronSection(x0, y0, z0, value0,
											   x1, y1, z1, value1,
											   x2, y2, z2, value2,
											   x3, y3, z3, value3);
		break;
	}
}

/*
 * Draw exact contour tetrahedron section
 */
static void drawExactContourTetrahedronSection(
	double x0, double y0, double z0, double value0,
	double x1, double y1, double z1, double value1,
	double x2, double y2, double z2, double value2,
	double x3, double y3, double z3, double value3) {
	if (contour_type == SMOOTH_CONTOUR_TYPE)
		AutoGL_DrawSectionSmoothContourTetrahedron(section_a,
				section_b,
				section_c,
				COORDS_SCALE_FACTOR * section_d,
				contour_field_min, contour_field_max,
				COORDS_SCALE_FACTOR * x0,
				COORDS_SCALE_FACTOR * y0,
				COORDS_SCALE_FACTOR * z0,
				value0,
				COORDS_SCALE_FACTOR * x1,
				COORDS_SCALE_FACTOR * y1,
				COORDS_SCALE_FACTOR * z1,
				value1,
				COORDS_SCALE_FACTOR * x2,
				COORDS_SCALE_FACTOR * y2,
				COORDS_SCALE_FACTOR * z2,
				value2,
				COORDS_SCALE_FACTOR * x3,
				COORDS_SCALE_FACTOR * y3,
				COORDS_SCALE_FACTOR * z3,
				value3);
	else if (contour_type == BAND_CONTOUR_TYPE)
		AutoGL_DrawSectionBandContourTetrahedron(section_a,
				section_b,
				section_c,
				COORDS_SCALE_FACTOR * section_d,
				contour_field_min, contour_field_max,
				num_contour_grades,
				COORDS_SCALE_FACTOR * x0,
				COORDS_SCALE_FACTOR * y0,
				COORDS_SCALE_FACTOR * z0,
				value0,
				COORDS_SCALE_FACTOR * x1,
				COORDS_SCALE_FACTOR * y1,
				COORDS_SCALE_FACTOR * z1,
				value1,
				COORDS_SCALE_FACTOR * x2,
				COORDS_SCALE_FACTOR * y2,
				COORDS_SCALE_FACTOR * z2,
				value2,
				COORDS_SCALE_FACTOR * x3,
				COORDS_SCALE_FACTOR * y3,
				COORDS_SCALE_FACTOR * z3,
				value3);
}

/*
 * Draw elemental contour tetrahedron section
 */
static void drawElementalContourTetrahedronSection(
	double x0, double y0, double z0, double value0,
	double x1, double y1, double z1, double value1,
	double x2, double y2, double z2, double value2,
	double x3, double y3, double z3, double value3) {
	double h0 = section_a * x0 + section_b * y0 + section_c * z0 + section_d;
	double h1 = section_a * x1 + section_b * y1 + section_c * z1 + section_d;
	double h2 = section_a * x2 + section_b * y2 + section_c * z2 + section_d;
	double h3 = section_a * x3 + section_b * y3 + section_c * z3 + section_d;

	int flag0 = 0;
	int flag1 = 0;
	int flag2 = 0;
	int flag3 = 0;

	switch (section_type) {
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
		flag0 = (h0 <= 0.0);
		flag1 = (h1 <= 0.0);
		flag2 = (h2 <= 0.0);
		flag3 = (h3 <= 0.0);
		break;
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		flag0 = (h0 >= 0.0);
		flag1 = (h1 >= 0.0);
		flag2 = (h2 >= 0.0);
		flag3 = (h3 >= 0.0);
		break;
	default:
		return;
	}

	if (flag0 && flag1 && flag2 && flag3)
		return;

	if (!flag0 && !flag1 && !flag2 && !flag3)
		return;

	if (flag1 && flag2 && flag3)
		drawContourTriangle(x1, y1, z1, value1,
							x2, y2, z2, value2,
							x3, y3, z3, value3);

	if (flag0 && flag3 && flag2)
		drawContourTriangle(x0, y0, z0, value0,
							x3, y3, z3, value3,
							x2, y2, z2, value2);

	if (flag0 && flag1 && flag3)
		drawContourTriangle(x0, y0, z0, value0,
							x1, y1, z1, value1,
							x3, y3, z3, value3);

	if (flag0 && flag2 && flag1)
		drawContourTriangle(x0, y0, z0, value0,
							x2, y2, z2, value2,
							x1, y1, z1, value1);
}

/*
 * Draw contour hexahedron section
 */
void drawContourHexahedronSection(double x0, double y0, double z0, double value0,
								  double x1, double y1, double z1, double value1,
								  double x2, double y2, double z2, double value2,
								  double x3, double y3, double z3, double value3,
								  double x4, double y4, double z4, double value4,
								  double x5, double y5, double z5, double value5,
								  double x6, double y6, double z6, double value6,
								  double x7, double y7, double z7, double value7) {
	switch (section_type) {
	case EQUAL_SECTION_TYPE:
	case LESS_EQUAL_SECTION_TYPE:
	case GREATER_EQUAL_SECTION_TYPE:
		drawExactContourHexahedronSection(x0, y0, z0, value0,
										  x1, y1, z1, value1,
										  x2, y2, z2, value2,
										  x3, y3, z3, value3,
										  x4, y4, z4, value4,
										  x5, y5, z5, value5,
										  x6, y6, z6, value6,
										  x7, y7, z7, value7);
		break;
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		drawElementalContourHexahedronSection(x0, y0, z0, value0,
											  x1, y1, z1, value1,
											  x2, y2, z2, value2,
											  x3, y3, z3, value3,
											  x4, y4, z4, value4,
											  x5, y5, z5, value5,
											  x6, y6, z6, value6,
											  x7, y7, z7, value7);
		break;
	}
}

/*
 * Draw exact contour hexahedron section
 */
static void drawExactContourHexahedronSection(
	double x0, double y0, double z0, double value0,
	double x1, double y1, double z1, double value1,
	double x2, double y2, double z2, double value2,
	double x3, double y3, double z3, double value3,
	double x4, double y4, double z4, double value4,
	double x5, double y5, double z5, double value5,
	double x6, double y6, double z6, double value6,
	double x7, double y7, double z7, double value7) {

	if (contour_type == SMOOTH_CONTOUR_TYPE)
		AutoGL_DrawSectionSmoothContourHexahedron(section_a,
				section_b,
				section_c,
				COORDS_SCALE_FACTOR * section_d,
				contour_field_min, contour_field_max,
				COORDS_SCALE_FACTOR * x0,
				COORDS_SCALE_FACTOR * y0,
				COORDS_SCALE_FACTOR * z0,
				value0,
				COORDS_SCALE_FACTOR * x1,
				COORDS_SCALE_FACTOR * y1,
				COORDS_SCALE_FACTOR * z1,
				value1,
				COORDS_SCALE_FACTOR * x2,
				COORDS_SCALE_FACTOR * y2,
				COORDS_SCALE_FACTOR * z2,
				value2,
				COORDS_SCALE_FACTOR * x3,
				COORDS_SCALE_FACTOR * y3,
				COORDS_SCALE_FACTOR * z3,
				value3,
				COORDS_SCALE_FACTOR * x4,
				COORDS_SCALE_FACTOR * y4,
				COORDS_SCALE_FACTOR * z4,
				value4,
				COORDS_SCALE_FACTOR * x5,
				COORDS_SCALE_FACTOR * y5,
				COORDS_SCALE_FACTOR * z5,
				value5,
				COORDS_SCALE_FACTOR * x6,
				COORDS_SCALE_FACTOR * y6,
				COORDS_SCALE_FACTOR * z6,
				value6,
				COORDS_SCALE_FACTOR * x7,
				COORDS_SCALE_FACTOR * y7,
				COORDS_SCALE_FACTOR * z7,
				value7);
	else if (contour_type == BAND_CONTOUR_TYPE)
		AutoGL_DrawSectionBandContourHexahedron(section_a,
												section_b,
												section_c,
												COORDS_SCALE_FACTOR * section_d,
												contour_field_min, contour_field_max,
												num_contour_grades,
												COORDS_SCALE_FACTOR * x0,
												COORDS_SCALE_FACTOR * y0,
												COORDS_SCALE_FACTOR * z0,
												value0,
												COORDS_SCALE_FACTOR * x1,
												COORDS_SCALE_FACTOR * y1,
												COORDS_SCALE_FACTOR * z1,
												value1,
												COORDS_SCALE_FACTOR * x2,
												COORDS_SCALE_FACTOR * y2,
												COORDS_SCALE_FACTOR * z2,
												value2,
												COORDS_SCALE_FACTOR * x3,
												COORDS_SCALE_FACTOR * y3,
												COORDS_SCALE_FACTOR * z3,
												value3,
												COORDS_SCALE_FACTOR * x4,
												COORDS_SCALE_FACTOR * y4,
												COORDS_SCALE_FACTOR * z4,
												value4,
												COORDS_SCALE_FACTOR * x5,
												COORDS_SCALE_FACTOR * y5,
												COORDS_SCALE_FACTOR * z5,
												value5,
												COORDS_SCALE_FACTOR * x6,
												COORDS_SCALE_FACTOR * y6,
												COORDS_SCALE_FACTOR * z6,
												value6,
												COORDS_SCALE_FACTOR * x7,
												COORDS_SCALE_FACTOR * y7,
												COORDS_SCALE_FACTOR * z7,
												value7);
}

/*
 * Draw elemental contour hexahedron section
 */
static void drawElementalContourHexahedronSection(
	double x0, double y0, double z0, double value0,
	double x1, double y1, double z1, double value1,
	double x2, double y2, double z2, double value2,
	double x3, double y3, double z3, double value3,
	double x4, double y4, double z4, double value4,
	double x5, double y5, double z5, double value5,
	double x6, double y6, double z6, double value6,
	double x7, double y7, double z7, double value7) {
	double h0 = section_a * x0 + section_b * y0 + section_c * z0 + section_d;
	double h1 = section_a * x1 + section_b * y1 + section_c * z1 + section_d;
	double h2 = section_a * x2 + section_b * y2 + section_c * z2 + section_d;
	double h3 = section_a * x3 + section_b * y3 + section_c * z3 + section_d;
	double h4 = section_a * x4 + section_b * y4 + section_c * z4 + section_d;
	double h5 = section_a * x5 + section_b * y5 + section_c * z5 + section_d;
	double h6 = section_a * x6 + section_b * y6 + section_c * z6 + section_d;
	double h7 = section_a * x7 + section_b * y7 + section_c * z7 + section_d;

	int flag0 = 0;
	int flag1 = 0;
	int flag2 = 0;
	int flag3 = 0;
	int flag4 = 0;
	int flag5 = 0;
	int flag6 = 0;
	int flag7 = 0;

	switch (section_type) {
	case ELEMENTAL_LESS_EQUAL_SECTION_TYPE:
		flag0 = (h0 <= 0.0);
		flag1 = (h1 <= 0.0);
		flag2 = (h2 <= 0.0);
		flag3 = (h3 <= 0.0);
		flag4 = (h4 <= 0.0);
		flag5 = (h5 <= 0.0);
		flag6 = (h6 <= 0.0);
		flag7 = (h7 <= 0.0);
		break;
	case ELEMENTAL_GREATER_EQUAL_SECTION_TYPE:
		flag0 = (h0 >= 0.0);
		flag1 = (h1 >= 0.0);
		flag2 = (h2 >= 0.0);
		flag3 = (h3 >= 0.0);
		flag4 = (h4 >= 0.0);
		flag5 = (h5 >= 0.0);
		flag6 = (h6 >= 0.0);
		flag7 = (h7 >= 0.0);
		break;
	default:
		return;
	}

	if (flag0 && flag1 && flag2 && flag3 && flag4 && flag5 && flag6 && flag7)
		return;

	if (!flag0 && !flag1 && !flag2 && !flag3 && !flag4 && !flag5 && !flag6 && !flag7)
		return;

	if (flag0 && flag3 && flag2 && flag1)
		drawContourQuadrangle(x0, y0, z0, value0,
							  x3, y3, z3, value3,
							  x2, y2, z2, value2,
							  x1, y1, z1, value1);

	if (flag4 && flag5 && flag6 && flag7)
		drawContourQuadrangle(x4, y4, z4, value4,
							  x5, y5, z5, value5,
							  x6, y6, z6, value6,
							  x7, y7, z7, value7);

	if (flag0 && flag1 && flag5 && flag4)
		drawContourQuadrangle(x0, y0, z0, value0,
							  x1, y1, z1, value1,
							  x5, y5, z5, value5,
							  x4, y4, z4, value4);

	if (flag1 && flag2 && flag6 && flag5)
		drawContourQuadrangle(x1, y1, z1, value1,
							  x2, y2, z2, value2,
							  x6, y6, z6, value6,
							  x5, y5, z5, value5);

	if (flag2 && flag3 && flag7 && flag6)
		drawContourQuadrangle(x2, y2, z2, value2,
							  x3, y3, z3, value3,
							  x7, y7, z7, value7,
							  x6, y6, z6, value6);

	if (flag0 && flag4 && flag7 && flag3)
		drawContourQuadrangle(x0, y0, z0, value0,
							  x4, y4, z4, value4,
							  x7, y7, z7, value7,
							  x3, y3, z3, value3);
}

/*
 * Draw point
 */
void drawPoint(double x, double y, double z) {
	switch (point_mark) {
	case CROSS_POINT_MARK:
		AutoGL_DrawCrossMark3D(COORDS_SCALE_FACTOR * x,
							   COORDS_SCALE_FACTOR * y,
							   COORDS_SCALE_FACTOR * z,
							   point_size);
		break;
	case DIAMOND_POINT_MARK:
		AutoGL_DrawDiamondMark3D(COORDS_SCALE_FACTOR * x,
								 COORDS_SCALE_FACTOR * y,
								 COORDS_SCALE_FACTOR * z,
								 point_size, 1);
		break;
	case CUBE_POINT_MARK:
		AutoGL_DrawCubeMark3D(COORDS_SCALE_FACTOR * x,
							  COORDS_SCALE_FACTOR * y,
							  COORDS_SCALE_FACTOR * z,
							  point_size, 1);
		break;
	}
}

/*
 * Draw mark
 */
void drawMark(double x, double y, double z) {
	AutoGL_DrawDiamondMark3D(COORDS_SCALE_FACTOR * x,
							 COORDS_SCALE_FACTOR * y,
							 COORDS_SCALE_FACTOR * z,
							 10, 1);
}
