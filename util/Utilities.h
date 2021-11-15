#ifndef _Utilities_h
#define _Utilities_h

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
//#include <CLHEP/Matrix/Vector.h>
//#include <CLHEP/Matrix/Matrix.h>
//#include <CLHEP/Matrix/SymMatrix.h>

/** \file Utilities.h
 *
 * \brief Various preprocessor macros, variables, functions, and procedures
 * that misfit elsewhere.
 *
 */

// global Variable to log errors (1), warnings (2), info (3), debug(4,5,...)
extern int gLogLevel;

#ifndef NDEBUG
#define LOG(level, message) if (gLogLevel >= level) { switch (level) { \
case 1: std::cerr << "ERROR: " << message << std::endl; break; \
case 2: std::cerr << "WARNING: " << message << std::endl; break; \
case 3: std::cout << "INFO: " << message << std::endl; break; \
default: std::cout << "DEBUG: " << message << std::endl; } }
#else
#define LOG(level, message) ;
#endif

#define ERROR(message) LOG(1, message);
#define WARNING(message) LOG(2, message);
#define INFO(message) LOG(3, message);

#define MATRIX_SIZE(mname, matrix) std::cout << "Matrix " << mname << " has " << matrix.num_row() << " rows and " << matrix.num_col() << " columns." << std::endl;
#define PMATRIX(mname, matrix) if (gLogLevel > 2) { cout << mname; matrix.Print(); }
#define MATRIX(matrix) PMATRIX(#matrix, matrix)

// throw an exception that tells me where the exception happened
#define THROW(errmsg) throw (std::string( __PRETTY_FUNCTION__ )+std::string(" (file: ")+std::string( __FILE__ )+std::string(", line: ")+std::string( Form("%d", __LINE__) )+std::string(") ")+std::string(errmsg));
#define CATCH catch (std::string message) { cerr << "EXCEPTION in " << message << std::endl; }
  
//void setXYZ(CLHEP::HepVector & v, double x, double y, double z);
bool equal(const double val1, const double val2, const double precision);

/** Get a submatrix from a large matrix */
//void GetSubMatrix(const CLHEP::HepMatrix & m, int startRow, int startColumn, int endRow, int endColumn, CLHEP::HepMatrix & sub);
//void GetSubMatrix(const CLHEP::HepSymMatrix & m, int startRow, int startColumn, int endRow, int endColumn, CLHEP::HepMatrix & sub);
/** Set a submatrix in a large matrix */
//void SetSubMatrix(CLHEP::HepMatrix & m, int startRow, int startColumn, const CLHEP::HepMatrix & sub);
//void SetSubMatrix(CLHEP::HepSymMatrix & m, int startRow, int startColumn, const CLHEP::HepMatrix & sub);

/** Get/set a subvector from a large vector */
//CLHEP::HepVector GetSubVector(const CLHEP::HepVector & v, int start, int end);
//void SetSubVector(CLHEP::HepVector & v, int start, const CLHEP::HepVector & sub);

//////////////////////////////////////////////////////////////////////
// Conversion of matrices

//void CLHEPtoROOT(const CLHEP::HepMatrix & oldM, TMatrixD * newM);
//void CLHEPtoROOT(const CLHEP::HepSymMatrix & oldM, TMatrixDSym * newM);
//void CLHEPtoROOT(const CLHEP::HepVector & oldV, TVectorD * newV);

//////////////////////////////////////////////////////////////////////
// Find objects from ROOT files

TObject * get_object(const char * objectname, const char * filename);
TTree * get_tree(const char * treename, const char * filename);
TTree * init_align_tree(const char * filename);

#endif
