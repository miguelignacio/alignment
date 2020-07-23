#include "KalmanFilterTracking.h"

// C++ includes
#include <iostream>

using namespace std;

// ROOT includes
#include "TMath.h"

// local includes
#include "Constants.h"
#include "Detector.h"
#include "Utilities.h"
#include "ThreeDModel.h"
#include "Visualization.h"

using namespace CLHEP;

/// Global Visualization object
extern Visualization gVis;

/** \file KalmanFilterTracking.cc
 *
 *  \brief Definition of KalmanFilterTracking class.
 */

KalmanFilterTracking::KalmanFilterTracking() : P(2,4,0)
{
  //////////////////////////////
  // setup projection matrices
  // projection from track local x and y to the 2d position
  P[0][2] = 1.;
  P[1][3] = 1.;  
}

KalmanFilterTracking::~KalmanFilterTracking()
{
}

/** \brief Jacobian matrix, derivatives of measurement equation to track
 *  parameters.
 *
 * The matrix is evaluated at the position p0, where p0 are the track
 * parameters.  p = (tx, ty, x, y)
 */
const HepMatrix & KalmanFilterTracking::Jacobian_Track(const HepVector & p0) const
{
  return P;
}

/** \brief A track is seeded from the first and last hit.
 * 
 * The seed is created from the first and from the last hit in the
 * hitvector. This already gives a good estimate of the track
 * parameters. However, the track is created with large initial covariance
 * matrix in order to keep the bias from the seed small.
 * 
 * When less than two hits are found on the track, CreateSeedTrack() is called
 * without argument.
 *
 */
Track KalmanFilterTracking::CreateSeedTrack(const std::vector<Hit *> & hitvector)
{
  assert(gDetector != 0);
  unsigned int size = hitvector.size();
  if (size < 2)
    return CreateSeedTrack();
  const Hit * firsthit = *hitvector.begin();
  const Hit * lasthit = *hitvector.rbegin();
  assert(firsthit != 0 && lasthit != 0);
  Track track(firsthit->GlobalMeasuredPosition(), lasthit->GlobalMeasuredPosition(), 4.);
  LOG(5, "Seed Track parameters" << track.GetGlobalParameters());
  /// assume errors for initial covariance matrix
  HepSymMatrix covariance(4, 0);
  for (int i = 0; i < 2; i++) {
    // there is a tradeoff between choosing large values (which makes the
    // tracking numerically unstable) and choosing small values (which biases
    // the track from the seed). I found 1E4 working for this model, but this
    // might need to get tuned.
    covariance[i][i] = 1E1;
    covariance[i+2][i+2] = 1E4;
  }
  track.SetGlobalCovariance(covariance);
  return track;
}

/** \brief A track seed is created in the origin going parallel to Z with huge
 * errors. */
Track KalmanFilterTracking::CreateSeedTrack()
{
  assert(gDetector != 0);

  /// Assume track originating in center and going along Z with very large error
  HepVector origin;
  double width  = gDetector->GetWidth();
  double length = gDetector->GetLength();
  double height = gDetector->GetHeight();
  setXYZ(origin, width/2., length/2., 0);
  HepVector destination;
  setXYZ(destination, width/2., length/2., height);
  Track track(origin, destination, 4.);
  HepSymMatrix covariance(4, 0);
  for (int i = 0; i < 2; i++) {
    // there is a tradeoff between choosing large values (which makes the
    // tracking numerically unstable) and choosing small values (which biases
    // the track from the seed). I found 1E4 working for this model, but this
    // might need to get tuned.
    covariance[i][i] = 1E1;
    covariance[i+2][i+2] = 1E4;
  }
  track.SetGlobalCovariance(covariance);
  return track;
}


/** \brief Process given hits and return a track built with Kalman Filter. 
 *
 * First, a track seed is built from the hit.
 *
 * Then, the track is updated by the Kalman Filter with each hit. The method
 * returns the final track after all hits are processed, with updated track
 * parameters and updated covariance matrix.
 *
 * In this method, no alignment is performed. The track is created taking into
 * account only the nominal position of the detectors, thus the alignment is
 * assumed not to be existing.
 *
 */
Track KalmanFilterTracking::MakeRecoTrack(vector<Hit *> & hitvector)
{
  Track track = CreateSeedTrack(hitvector);
  /// Now process hits and estimate track parameters
  int i = 0;
  for (vector<Hit *>::const_iterator it = hitvector.begin(); it != hitvector.end(); it++, i++) {
    const Hit & hit = *(*it);
    const Det & det = hit.GetDet();
    // LOG(7, "KFA::MakeReferenceTrack(): hit in layer " << det.GetLayer() << " rod " << det.GetRod() << " module " << det.GetModule());
    const ReferenceFrame frame = det.GetNominal();

    /// previous track state
    LOG(7, "Initial track state (global) " << track.GetGlobalParameters());
    HepVector p0 = track.GetLocalParameters(frame);
    LOG(7, "Initial track state (local)" << p0);
    /// actual measurement (it is only 2-dimensional in x and y)
    HepVector m  = GetSubVector(hit.GetPosition(), 0, 1);
    LOG(7, "actual measurement" << m);
    /// predicted measurement given previous track state
    HepVector f0 = P*p0;
    LOG(7, "Predicted measurement" << f0);

    // split into derivatives to track parameters
    HepMatrix H = Jacobian_Track(p0);
    LOG(7, "Track derivatives" << H);

    HepSymMatrix V = hit.GetCovariance();
    LOG(7, "hit covariance matrix " << V);
    HepSymMatrix C0 = track.GetLocalCovariance(frame);
    LOG(7, "track covariance matrix C0 " << C0);
    // compute the weight matrix
    int ierr;
    HepSymMatrix W = (V + C0.similarity(H)).inverse(ierr);
    if (ierr != 0)
      THROW("matrix inversion failed");

    LOG(7, "Weight matrix " << W);

    // residual estimate vs. measurement
    HepVector r = m - f0;
    LOG(7, "residual " << r);

    // \todo here one should add protection against adding very bad hits to
    // the track, based on the current track and hit covariance matrix and the
    // residual. Problems can arise due to two reasons: Either the track seed
    // is bad (i.e. points in a wrong direction with small errors), or the
    // hits are from detector noise or other tracks.

    // Kalman weight matrix K
    HepMatrix K = C0 * H.T() * W;
    LOG(7, "Kalman weight matrix track " << K);
    // track parameter update
    HepVector p1 = p0 +  K * r;
    LOG(7, "new track parameters " << p1);
    HepSymMatrix C1 = C0 - (W.similarityT(H)).similarity(C0);
    LOG(7, "new track covariance " << C1);

    // set new track parameters
    track.SetLocalParameters(frame, p1);
    track.SetLocalCovariance(frame, C1);

    LOG(7, "new track global parameters" << track.GetGlobalParameters());
    LOG(7, "Track parameter errors" << endl
	<< TMath::Sqrt(track.GetGlobalCovariance()[0][0]) << endl
	<< TMath::Sqrt(track.GetGlobalCovariance()[1][1]) << endl
	<< TMath::Sqrt(track.GetGlobalCovariance()[2][2]) << endl
	<< TMath::Sqrt(track.GetGlobalCovariance()[3][3]) << endl);
    // // show track, but exclude first measurement (because tangent undefined)
    // if (i > 0)
    //   gVis.AddTrack(track, i+1);

    // save hit for future access
    // We have to duplicate hit since old hits are already owned.
    track.GetHits().push_back(new Hit(*(*it)));
  }
  // now that the track has been reconstructed, we will compute the track
  // chi^2.  Therefore it is necessary to loop once again over the hits to
  // evaluate the track again at each point
  HepVector chi2(1, 0);
  for (vector<Hit *>::const_iterator it = hitvector.begin(); it != hitvector.end(); it++) {
    const Hit & hit = *(*it);
    const Det & det = hit.GetDet();
    Hit * trackPrediction = track.NominalIntersection(det);
    if (trackPrediction == 0) {
      // \todo It would be better not to throw away the full track but rather
      // only the hit. In some configuration, this happens for about 1 percent
      // of the tracks. For this simulation efficiency does not matter,
      // therefore I do not care right now.
      THROW("This track has no hit on the given Det.");
      continue;
    }
    HepSymMatrix trackErrors = track.GetLocalCovariance(det.GetNominal()).similarity(P);
    HepVector residual = (hit.GetPosition()-trackPrediction->GetPosition()).sub(1,2);
    HepSymMatrix errorMatrix = hit.GetCovariance();
    assert((errorMatrix[0][0] > 0) && (errorMatrix[1][1] > 0));
    int ierr;
    chi2 += residual.T()*(errorMatrix.inverse(ierr))*residual;
    if (ierr != 0)
      THROW("matrix inversion failed");
  }
  track.SetChi2Ndof(chi2[0], hitvector.size()*2-4);
//   gVis.AddTrack(track);
  return track;
}
