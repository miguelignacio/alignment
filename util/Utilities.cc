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

/** \file Utilities.h
 *
 * \brief Various preprocessor macros, variables, functions, and procedures
 * that misfit elsewhere.
 *
 */

// global logging level
int gLogLevel = 3;


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

