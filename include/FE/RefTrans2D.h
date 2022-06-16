// =======================================================================
// @(#)RefTrans2D.h        1.4 04/13/00
//
// Class:      TRefTrans2D
//
// Purpose:    reference transformations for 2D geometric objects
//
// Author:     Gunar Matthies
//
// History:    08.07.97 start implementation
//
// =======================================================================

#ifndef __REFTRANS2D__
#define __REFTRANS2D__

#include <Constants.h>
#include <Enumerations.h>
#include <BaseCell.h>

/** reference transformations for 2D geometric objects */
class TRefTrans2D
{
protected:
	TBaseCell *Cell;

public:
	/** constuctor */
	TRefTrans2D(){};

	/** transfer form reference element to original element */
	void GetOrigFromRef(double eta, double xi, double &x, double &y);

	/** transfer form reference element to original element */
	void GetOrigFromRef(double *ref, double *orig);

	/** transfer from original element to reference element */
	virtual void GetRefFromOrig(double x, double y, double &eta, double &xi)
	{
		cout << "[ERROR - Thivin] BaseClass polymorphism Error.  File : RefTrans2D.h" << endl;
		cout << " Unable to call the corresponding child class for reference transformation. This is an virtual fucntion. Please make sure the corresponding child class has same function with the definition" << endl;
		exit(0);
	}

	/** transfer from original element to reference element */
	void GetRefFromOrig(double *orig, double *ref);

	/** calculate functions and derivatives from reference element
        to original element */
	void GetOrigValues(TBaseCell *cell);

	/** Get original Values from the Reference element for given values of shape functions and derivatives at given points **/
	// Written by Thivin
	virtual void GetOrigValues(double xi, double eta,
							   int N_BaseFunct,
							   double *uref, double *uxiref, double *uetaref,
							   double *uorig, double *uxorig, double *uyorig,
							   int _BaseVectDim)
	{
		cout << "[ERROR]-Thivin File: RefTrans2D.h " <<endl;
		cout << " \"GetOrigValues\" Not implemented for the given Reftranstype. Either check the implementation in the respective class" <<endl;
		exit(0);
	}

	/** set original element to cell */
	virtual void SetCell(TBaseCell *cell)
	{
		cout << " Setting cell from the Base class -- Not expected behaviour, Check implementation  " << endl;
		Cell = cell;
	}

	static RefTrans2D FindRefTrans2D(int N_LocalUsedElements, FE2D *LocalUsedElements);

	/** return outer normal vector */
	virtual void GetOuterNormal(int j, double zeta,
								double &n1, double &n2) = 0;

	/** return tangent */
	virtual void GetTangent(int j, double zeta,
							double &t1, double &t2) = 0;

	/** return volume of cell according to reference transformation */
	double GetVolume();

	// piola map, needed for vector values basis functions such as
	// Raviart-Thomas (RT) or Brezzi-Douglas-Marini (BDM)
	virtual void PiolaMapOrigFromRef(int N_Functs, double *refD00, double *origD00)
	{
		cout << " Piola Map not defined for this element " << endl;
	};

	// piola map for derivatives
	virtual void PiolaMapOrigFromRef(int N_Functs, double *refD10,
									 double *refD01, double *origD10,
									 double *origD01)
	{
		cout << " Piola Map not defined for this element " << endl;
	};

	virtual ~TRefTrans2D()
	{
	}
};

#endif
