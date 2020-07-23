// C++ includes
#include <iostream>
#include <assert.h>

using namespace std;

// CLHEP includes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>

using namespace CLHEP;

// local includes
#include "Detector.h"
#include "Hit.h"
#include "Utilities.h"

/** 
 * \file Hit.cc 
 *
 * \brief Definition of class Hit. 
 * 
 * Internally, the hit coordinates are saved as local coordinates. Therefore,
 * the z hit coordinate is always zero.
 *
 */

Hit::Hit(const Det & det, HepVector position) : fDet(&det), fPosition(position)
{
}

Hit::Hit(const Det & det) : fDet(&det)
{
  // set position
  setXYZ(fPosition, 0, 0, 0);
  // check that the given position is inside detector bounds
  assert(fDet->InsideBounds(fPosition));
}

/** \brief Return the global measured position of the hit. 
 *
 * I.e. compute the global position with true misaligned position. 
 */
HepVector Hit::GlobalMeasuredPosition() const
{
  return fDet->GetMisalign().ToGlobal(fPosition);
}

/** Return covariance matrix of hit */
HepSymMatrix Hit::GetCovariance() const
{
  /// Return covariance matrix
  HepSymMatrix C(2, 0);
  double xRes = fDet->GetResolutionX();
  double yRes = fDet->GetResolutionY();
  C[0][0] = xRes*xRes;
  C[1][1] = yRes*yRes;

  return C;
}

HepVector Hit::GlobalTruePosition() const
{
  return fDet->GetMisalign().ToGlobal(fPosition);
}
