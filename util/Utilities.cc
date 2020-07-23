#include "Utilities.h"

/** \file Utilities.cc Contains often used functions for vectors,  matrices and other. */

// C++ includes
#include <assert.h>
#include <iostream>

// ROOT includes
#include <TMath.h>
#include <TFile.h>
#include <TROOT.h>

using namespace std;
using namespace CLHEP;

/** \file Utilities.h
 *
 * \brief Various preprocessor macros, variables, functions, and procedures
 * that misfit elsewhere.
 *
 */

// global logging level
int gLogLevel = 3;

/** Set X, Y and Z coordinate of a HepVector in one go. */
void setXYZ(HepVector & v, double x, double y, double z)
{
  if (v.num_row() != 3) {
//     cout << "Resizing vector from size " << v.num_row() << " to size 3" << endl;
    HepVector t(3);
    v = t;
  }
  v[0] = x;
  v[1] = y;
  v[2] = z;
}

/** Test if two values are equal within a certain precision. */
bool equal(const double val1, const double val2, const double precision)
{
  if (val1 == 0)
    return val2 < precision;
  if (val2 == 0)
    return val1 < precision;
//   cout << "val1= " << val1 << "   val2= " << val2 << "  equality= " << TMath::Abs(2*(val1-val2)/(val1+val2)) << endl;
  return TMath::Abs(2*(val1-val2)/(val1+val2)) < precision;
}

/** Get a submatrix of matrix m. Counting starts at 0, and both start and end
 * rows/columns are included. The matrix sub needs to have the correct size,
 * already. */
void GetSubMatrix(const HepSymMatrix & m, int startRow, int startColumn, int endRow, int endColumn, HepMatrix & sub)
{
//   cout << "(test1, " << flush;
//   HepMatrix sub(endRow-startRow + 1, endColumn-startColumn+1);
  assert(sub.num_row() == endRow-startRow+1);
  assert(sub.num_col() == endColumn-startColumn+1);
  assert(endRow < m.num_row());
  assert(endColumn < m.num_col());
  for (int i = startRow; i <= endRow; i++) {
    for (int j = startColumn; j <= endColumn; j++) {
      sub(i-startRow+1,j-startColumn+1) = m(i+1,j+1);
    }
  }
//   cout << "test2)" << flush;
//   return sub;
}

/** Get a submatrix of matrix m. Counting starts at 0, and both start and end
 * rows/columns are included. The matrix sub needs to have the correct size,
 * already. */
void GetSubMatrix(const HepMatrix & m, int startRow, int startColumn, int endRow, int endColumn, HepMatrix & sub)
{
//   cout << "(test1, " << flush;
//   HepMatrix sub(endRow-startRow + 1, endColumn-startColumn+1);
  assert(sub.num_row() == endRow-startRow+1);
  assert(sub.num_col() == endColumn-startColumn+1);
  assert(endRow < m.num_row());
  assert(endColumn < m.num_col());
  for (int i = startRow; i <= endRow; i++) {
    for (int j = startColumn; j <= endColumn; j++) {
      sub(i-startRow+1,j-startColumn+1) = m(i+1,j+1);
    }
  }
//   cout << "test2)" << flush;
//   return sub;
}


/** Set a submatrix of the covariance matrix. Counting starts at 0. */
void SetSubMatrix(HepMatrix & m, int startRow, int startColumn, const HepMatrix & sub)
{
  assert(startRow >= 0);
  assert(startColumn >= 0);
  assert(startRow + sub.num_row() <= m.num_row());
  assert(startColumn + sub.num_col() <= m.num_col());
  for (int i = 1; i <= sub.num_row(); i++) {
    for (int j = 1; j <= sub.num_col(); j++) {
      m(i+startRow,j+startColumn) = sub(i,j);
    }
  }
}

/** Set a submatrix of the covariance matrix. Counting starts at 0. */
void SetSubMatrix(HepSymMatrix & m, int startRow, int startColumn, const HepMatrix & sub)
{
  assert(startRow >= 0);
  assert(startColumn >= 0);
  assert(startRow + sub.num_row() <= m.num_row());
  assert(startColumn + sub.num_col() <= m.num_col());
  for (int i = 1; i <= sub.num_row(); i++) {
    for (int j = 1; j <= sub.num_col(); j++) {
      m(i+startRow,j+startColumn) = sub(i,j);
    }
  }
}

HepVector GetSubVector(const HepVector & v, int start, int end)
{
  assert(start >= 0);
  assert(end < v.num_row());
  HepVector t(end-start+1);
  for (int i = start; i <= end; i++) {
    t[i-start] = v[i];
  }
  return t;
}

void SetSubVector(HepVector & v, int start, const HepVector & sub)
{
  assert(start >= 0);
  assert(start+sub.num_row() <= v.num_row());
  for (int i = 0; i < sub.num_row(); i++) {
    v[start+i] = sub[i];
  }
}

//////////////////////////////////////////////////////////////////////
// Conversion of matrices

/** \brief Convert a CLHEP HepMatrix to ROOT TMatrixD. */
void CLHEPtoROOT(const HepMatrix & oldM, TMatrixD * newM)
{
  newM->ResizeTo(oldM.num_row(), oldM.num_col());
  for (int i = 0; i < oldM.num_row(); i++) {
    for (int j = 0; j < oldM.num_col(); j++) {
      (*newM)[i][j] = oldM[i][j];
    }
  }
}

/** \brief Convert a CLHEP HepSymMatrix to ROOT TMatrixDSym. */
void CLHEPtoROOT(const HepSymMatrix & oldM, TMatrixDSym * newM)
{
//   cout << "YES!" << endl;
  newM->ResizeTo(oldM.num_row(), oldM.num_col());
  for (int i = 0; i < oldM.num_row(); i++) {
    for (int j = i; j < oldM.num_col(); j++) {
      (*newM)[i][j] = oldM[i][j];
    }
  }
}

/** \brief Convert a CLHEP HepVector to ROOT TVectorD. */
void CLHEPtoROOT(const HepVector & oldV, TVectorD * newV)
{
  newV->ResizeTo(oldV.num_row());
  for (int i = 0; i < oldV.num_row(); i++) {
      (*newV)[i] = oldV[i];
  }
}

// Utility routines
//==================

TObject * get_object(const char * objectname, const char * filename)
{
  // try to find out if file is already opened
  TFile * f = (TFile *) gROOT->GetListOfFiles()->FindObject(filename);
  if (f != 0) {
//     cout << "File " << fname << " was already open" << endl;
  }
  else {
    f = new TFile(filename, "READ");
    if (f == 0) {
      cout << "Error creating TFile object" << endl;
      return 0;
    }
    if (!f->IsOpen()) {
      cout << "Could not open file " << filename << endl;
      return 0;
    }
  }
  TObject * obj = (TTree *) f->Get(objectname);
  if (obj == 0) {
    cout << "Error reading object " << objectname << " from file " << filename << endl;
  }
  return obj;
}

TTree * get_tree(const char * treename, const char * filename)
{
  return (TTree *) get_object(treename, filename);
}

TTree * init_align_tree(const char * filename)
{
  TTree * t = get_tree("AlignTree", filename);
  return t;
}

