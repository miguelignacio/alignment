#ifndef _KalmanFilterTracking_h
#define _KalmanFilterTracking_h

// STL includes
#include <vector>

// CLHEP includes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>

// local includes
#include "Hit.h"
#include "Track.h"

/** \file KalmanFilterTracking.h
 *
 *  \brief Declaration of KalmanFilterTracking class.
 */

/** \class KalmanFilterTracking
 *  \brief Forward Kalman filter tracking algorithm, without smoothing step.
 *
 * This class allows straight track Kalman Filter tracking. This is done after
 * pattern recognition, i.e. the hits that constitute the track have to be
 * given to the algorithm.
 *
 */

class KalmanFilterTracking {
 protected:
  /// Create a seed track (from the given hits)
  Track CreateSeedTrack();
  Track CreateSeedTrack(const std::vector<Hit *> & hitvector);

 public:
  KalmanFilterTracking();
  virtual ~KalmanFilterTracking();

  /// jacobian matrix, derivatives of measurement equation to track parameters
  /// evaluated at initial track parameters p0
  const CLHEP::HepMatrix & Jacobian_Track(const CLHEP::HepVector & p0) const;

  /// reconstruct a track with KF from given hits
  Track MakeRecoTrack(std::vector<Hit *> & hitvector);

protected:
  /// projection matrix from track parameters to local measurement
  CLHEP::HepMatrix P;
};

#endif
