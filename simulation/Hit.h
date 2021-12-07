#ifndef _Hit_h
#define _Hit_h

#include "SimUtils.h"

// CLHEP includes
#include <CLHEP/Matrix/Vector.h>

#include "Det.h"

/** \file Hit.h 
 *  \brief Declaration of class Hit.
 */

/** \class Hit \brief A hit in a detector. Contains position, covariance
 *  matrix and link to Det. The hit position is given in local coordinates and
 *  can be retrieved in local or global coordinates.
 */

class Hit {
  /// Detector which owns the hit
  const Det * fDet;
  /// position of hit in local coordinates
  CLHEP::HepVector fPosition;
 public:
  Hit(const Det & det, CLHEP::HepVector position);
  Hit(const Det & det);
  class HitNotFoundException {};
  /// Detector which owns the hit
  void SetDet(const Det & det) { fDet = & det; };
  const Det & GetDet() const { return *fDet; };
  /// set position of hit in local coordinates
  void SetPosition(CLHEP::HepVector position) { fPosition = position; };
  CLHEP::HepVector GetPosition() const { return fPosition; };
  /// get position of hit as it would be measured
  CLHEP::HepVector GlobalMeasuredPosition() const;
  CLHEP::HepVector GlobalTruePosition() const;
  CLHEP::HepSymMatrix GetCovariance() const;
};

#endif
