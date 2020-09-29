// =======================================================================
// SSMUM.h
//
// Purpose:     shear slip mesh update method
//
// Authors:     Volker John, 2008/05/22
//
// =======================================================================

#ifndef __SSMUM__
#define __SSMUM__

/*******************************************************************************/
//
// WriteGridGnu
// writes a gnuplot output for the grid, just for debugging
//
/*******************************************************************************/

int WriteGridGnu(const char *name, TCollection *coll);


/*******************************************************************************/
//
// checks if the point (x,y) is in cell 
// works only for triangles
// input: cell and x,y coordinates of the vertices
// output 0 - no, 1 yes
// 
/*******************************************************************************/

int PointInCell(TBaseCell *cell, double x_coord, double y_coord);

/*******************************************************************************/
//
// SwapEdges
// routine for swapping the edges in the shear slip mesh update method in 2d
//
/*******************************************************************************/

void SwapEdges(const char *name, TCollection *coll,
	       TFEFunction2D *u1, TFEFunction2D *u2, 
	       double *tangential_values_ssl);

/*******************************************************************************/
//
// RotateGrid
// controls the grid rotation
//
/*******************************************************************************/

int RotateGrid(const char *name, TCollection *coll, double swap_rotation, double &angle,
	       double *uoldx, double *uoldy, double *uold1, double *uold2,
	       TFEFunction2D *u1, TFEFunction2D *u2, double *tangential_values_ssl);


/*******************************************************************************/
//
// VelocityInNewPositions
// computes the velicity values in the position after the rotation 
//
/*******************************************************************************/

int VelocityAtNewPositions(TCollection *coll, 
			   TFEFunction2D *u1, TFEFunction2D *u2,
			   double *values);

void FillNewVelocity(TCollection *coll, 
		     double *uoldx, double *uoldy, double *uold1, double *uold2,
		     TFEFunction2D *u1, TFEFunction2D *u2, double *tangential_values_ssl);

void MakeBubblesDivFree(TCollection *coll, 
			TFEFunction2D *u1, TFEFunction2D *u2, 
			TFEFunction2D *p, 
			TMatrix2D *matrixB1, TMatrix2D *matrixB2);

#endif
