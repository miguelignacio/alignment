#ifndef _SimUtil_h
#define _SimUtil_h

// C++ / STL includes
#include <string>
#include <iostream>
#include <assert.h>

// ROOT includes
#include "TString.h" // for char *Form(...)
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TObject.h"
#include "TTree.h"

// CLHEP includes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>

#include "Utilities.h"

/** \file SimUtil.h
 *
 * \brief Various preprocessor macros, variables, functions, and procedures
 * that misfit elsewhere.  All CLHEP macros, which are not used in the alignment algorithm,
 * have been moved to this module.
 *
 */

  
void setXYZ(CLHEP::HepVector & v, double x, double y, double z);

/** Get a submatrix from a large matrix */
void GetSubMatrix(const CLHEP::HepMatrix & m, int startRow, int startColumn, int endRow, int endColumn, CLHEP::HepMatrix & sub);
void GetSubMatrix(const CLHEP::HepSymMatrix & m, int startRow, int startColumn, int endRow, int endColumn, CLHEP::HepMatrix & sub);
/** Set a submatrix in a large matrix */
void SetSubMatrix(CLHEP::HepMatrix & m, int startRow, int startColumn, const CLHEP::HepMatrix & sub);
void SetSubMatrix(CLHEP::HepSymMatrix & m, int startRow, int startColumn, const CLHEP::HepMatrix & sub);

/** Get/set a subvector from a large vector */
CLHEP::HepVector GetSubVector(const CLHEP::HepVector & v, int start, int end);
void SetSubVector(CLHEP::HepVector & v, int start, const CLHEP::HepVector & sub);

//////////////////////////////////////////////////////////////////////
// Conversion of matrices

void CLHEPtoROOT(const CLHEP::HepMatrix & oldM, TMatrixD * newM);
void CLHEPtoROOT(const CLHEP::HepSymMatrix & oldM, TMatrixDSym * newM);
void CLHEPtoROOT(const CLHEP::HepVector & oldV, TVectorD * newV);


#endif
