#ifndef _KalmanFilterAlignmentInputProvider_h
#define _KalmanFilterAlignmentInputProvider_h

// CLHEP includes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include "SimUtils.h"

// local includes
#include <AlignEvent.h>
#include "Track.h"

/** \file KalmanFilterAlignmentInputProvider.h
 *
 *  \brief Declaration of KalmanFilterAlignmentInputProvider class.
 */

/** \class KalmanFilterAlignmentInputProvider
 *  \brief Take a track and produce the input for the Kalman Filter 
 *         Alignment Algorithm.
 *
 *  Once a reconstructed track is available, it can be used for
 *  alignment. This class prepares the input the Kalman Filter Alignment
 *  Algorithm will need to perform alignment. However, it is more general: The
 *  input can server for other alignment algorithms as well.
 *
 */

class KalmanFilterAlignmentInputProvider {
protected:
  const Track & fTrack; // the track from which information has to be given

  /// jacobian matrix, derivatives of measurement equation to alignment parameters, evaluated at p0 (and a0)
  const CLHEP::HepMatrix Jacobian_Alignment(const CLHEP::HepVector & p0) const;
  const CLHEP::HepMatrix Jacobian_Alignment(const CLHEP::HepVector & p0,
					    const CLHEP::HepVector & a0) const;

public:
  KalmanFilterAlignmentInputProvider(const Track & recoTrack);
  virtual ~KalmanFilterAlignmentInputProvider();

  /// compute all necessary values and fill Event structure
  void Fill(AlignEvent & event);
};

#endif
