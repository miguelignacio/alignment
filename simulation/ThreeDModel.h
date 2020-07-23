#ifndef _ThreeDModel_h
#define _ThreeDModel_h

// CLHEP includes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>

/** \file ThreeDModel.h
 *
 *  \brief Declaration of functions and procedures for 3D manipulation with matrices.
 */

// Definition of angles and rotations from Veikko Karimaki
void FillRotMatrixKarimaki(CLHEP::HepMatrix & m, double alpha, double beta, double gamma);
void GetAnglesKarimaki(const CLHEP::HepMatrix & m, double & alpha, double & beta, double & gamma);
void testAnglesKarimaki();

// Definition of angles  and rotations in polar coordinates
void FillRotMatrixPolar(CLHEP::HepMatrix & m, double theta, double phi);
// void GetAnglesPolar(const CLHEP::HepMatrix & m, double & theta, double & phi); // \todo implement.
void testAnglesPolar();

// test mixed angles
/* void testMixedAngles(int num); */

#endif 
