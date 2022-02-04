#ifndef DRAW_H
#define DRAW_H

/*
 * Draw line
 */
void drawLine(double x0, double y0, double z0,
			  double x1, double y1, double z1);

/*
 * Draw sectioned line
 */
void drawSectionedLine(double x0, double y0, double z0,
					   double x1, double y1, double z1);

/*
 * Draw tetrahedron section lines
 */
void drawTetrahedronSectionLines(double x0, double y0, double z0,
								 double x1, double y1, double z1,
								 double x2, double y2, double z2,
								 double x3, double y3, double z3);

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
								double x7, double y7, double z7);

/*
 * Draw triangle
 */
void drawTriangle(double x0, double y0, double z0,
				  double x1, double y1, double z1,
				  double x2, double y2, double z2);

/*
 * Draw quadrangle
 */
void drawQuadrangle(double x0, double y0, double z0,
					double x1, double y1, double z1,
					double x2, double y2, double z2,
					double x3, double y3, double z3);

/*
 * Draw sectioned triangle
 */
void drawSectionedTriangle(double x0, double y0, double z0,
						   double x1, double y1, double z1,
						   double x2, double y2, double z2);

/*
 * Draw sectioned quadrangle
 */
void drawSectionedQuadrangle(double x0, double y0, double z0,
							 double x1, double y1, double z1,
							 double x2, double y2, double z2,
							 double x3, double y3, double z3);

/*
 * Draw contour triangle
 */
void drawContourTriangle(double x0, double y0, double z0, double value0,
						 double x1, double y1, double z1, double value1,
						 double x2, double y2, double z2, double value2);

/*
 * Draw contour quadrangle
 */
void drawContourQuadrangle(double x0, double y0, double z0, double value0,
						   double x1, double y1, double z1, double value1,
						   double x2, double y2, double z2, double value2,
						   double x3, double y3, double z3, double value3);

/*
 * Draw sectioned contour triangle
 */
void drawSectionedContourTriangle(double x0, double y0, double z0, double value0,
								  double x1, double y1, double z1, double value1,
								  double x2, double y2, double z2, double value2);

/*
 * Draw sectioned contour quadrangle
 */
void drawSectionedContourQuadrangle(double x0, double y0, double z0, double value0,
									double x1, double y1, double z1, double value1,
									double x2, double y2, double z2, double value2,
									double x3, double y3, double z3, double value3);

/*
 * Draw contour tetrahedron section
 */
void drawContourTetrahedronSection(double x0, double y0, double z0, double value0,
								   double x1, double y1, double z1, double value1,
								   double x2, double y2, double z2, double value2,
								   double x3, double y3, double z3, double value3);

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
								  double x7, double y7, double z7, double value7);

/*
 * Draw point
 */
void drawPoint(double x, double y, double z);

/*
 * Draw mark
 */
void drawMark(double x, double y, double z);

#endif  /* DRAW_H */
