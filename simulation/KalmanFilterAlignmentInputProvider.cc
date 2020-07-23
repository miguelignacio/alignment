#include "KalmanFilterAlignmentInputProvider.h"

// C++ includes
#include <iostream>

using namespace std;

// ROOT includes
#include "TMath.h"

// local includes
#include "Constants.h"
#include "Utilities.h"
#include "KalmanFilterTracking.h"

using namespace CLHEP;

/** \file KalmanFilterAlignmentInputProvider.cc
 *
 *  \brief Definition of KalmanFilterAlignmentInputProvider class.
 */

/** \brief Constructor.
 */
KalmanFilterAlignmentInputProvider::KalmanFilterAlignmentInputProvider(const Track & recoTrack) :
  fTrack(recoTrack)
{
}

/** \brief Destructor.
 */
KalmanFilterAlignmentInputProvider::~KalmanFilterAlignmentInputProvider()
{
  // nothing to do!
}

/// compute all necessary matrices
void KalmanFilterAlignmentInputProvider::Fill(AlignEvent & event)
{
  // number of hits on track
  const int nHits = fTrack.GetHits().size();
  // number of track parameters
  const int nTrackParam = fTrack.GetGlobalParameters().num_row();
  // instantiate a tracking class (needed for derivatives)
  KalmanFilterTracking tracking;

  //////////////////////////////////////////////////////////////////////
  // create and initialize vectors and matrices to be used

  // construct ID array.
  int * ids = new int[nHits];

  // matrix of derivatives to track parameters
  HepMatrix TrackDerivatives(2*nHits, nTrackParam, 0);

  // matrix of derivatives to alignment parameters
  HepMatrix AlignmentDerivatives(2*nHits, nAlignmentParam*nHits, 0);

  // predicted measurement
  HepVector TrackPrediction(2*nHits, 0);

  // actual measurement
  HepVector Measurement(2*nHits, 0);

  // measurement covariance matrix
  HepSymMatrix MeasuredCovariance(2*nHits, 0);

  //////////////////////////////////////////////////////////////////////
  // fill local vectors and matrices with values

  // for all other values, a loop over detectors is necessary
  for (int i = 0; i < nHits; i++) {
    const Hit & hit = *(fTrack.GetHits()[i]);
    const Det & det = hit.GetDet();

    // fill ID
    ids[i] = det.GetID();

    // get prediction in nominal frame
    const ReferenceFrame frame = det.GetNominal();

    // initial track state
    HepVector p0 = fTrack.GetLocalParameters(frame);
//     cout << "Initial track state (local)" << p0;

    // actual measurement
    SetSubVector(Measurement, i*2, GetSubVector(hit.GetPosition(), 0, 1));

    // Jacobian matrix, derivatives of measurement equation to local track parameters
    const HepMatrix Jpl = tracking.Jacobian_Track(p0);

    // predicted measurement
    HepVector tp = Jpl*p0;
    SetSubVector(TrackPrediction, i*2, tp);

    // need to convert Jacobian for track parameters (local) into the global
    // jacobian.  this is done by multiplying with the local->global
    // transformation matrix
    const HepMatrix Jpg = Jpl * fTrack.Jacobian(frame);
//     cout << "global track derivative matrix" << Jpg;
    SetSubMatrix(TrackDerivatives, 2*i, 0, Jpg);

    // Jacobian matrix, derivatives of measurement equation to alignment parameters
    const HepMatrix Ja = Jacobian_Alignment(p0);
    // fill in derivatives to alignment parameters
    SetSubMatrix(AlignmentDerivatives, 2*i, nAlignmentParam*i, Ja);

    // measurement covariance matrix
    SetSubMatrix(MeasuredCovariance, i*2, i*2, hit.GetCovariance());
  }
  // MATRIX(TrackPrediction);
  // MATRIX(Measurement);
  // MATRIX(MeasuredCovariance);
  // MATRIX(TrackDerivatives);
  // MATRIX(AlignmentDerivatives);

  //////////////////////////////////////////////////////////////////////
  // Fill event structure
  event.SetChi2(fTrack.GetChi2());
  event.SetNdof(fTrack.GetNdof());
  // let TArrayI take ownership of the ids...
  event.GetIndex()->Adopt(nHits, ids);
  // save all matrices in tree
  CLHEPtoROOT(Measurement, event.GetMeasurements());
  CLHEPtoROOT(MeasuredCovariance, event.GetMeasuredCovariance());
  CLHEPtoROOT(TrackDerivatives, event.GetTrackDerivatives());
  CLHEPtoROOT(AlignmentDerivatives, event.GetAlignmentDerivatives());
  CLHEPtoROOT(TrackPrediction, event.GetTrackPrediction());
}

/** \brief Jacobian matrix, derivatives of measurement equation to alignment
 *  parameters.
 *
 *  The matrix is evaluated at the position p0, a0, where p0 and a0 are the
 *  track and alignment parameters, respectively.  p = (tx, ty, x, y), a =
 *  (dx, dy, dz, dalpha, dbeta, dgamma)
 */
const HepMatrix KalmanFilterAlignmentInputProvider::Jacobian_Alignment(const HepVector & p0, 
								       const HepVector & a0) const
{
  // the jacobian matrix
  HepMatrix J(2, 6);

  // prepare parameters
  const double tx = p0[0];
  const double ty = p0[1];
  const double x  = p0[2];
  const double y  = p0[3];
  const double dx = a0[0];
  const double dy = a0[1];
  const double dz = a0[2];
  const double sdalpha = TMath::Sin(a0[3]);
  const double cdalpha = TMath::Cos(a0[3]);
  const double sdbeta  = TMath::Sin(a0[4]);
  const double cdbeta  = TMath::Cos(a0[4]);
  const double sdgamma = TMath::Sin(a0[5]);
  const double cdgamma = TMath::Cos(a0[5]);

  const double norm = (cdalpha*cdbeta + ty*cdbeta*sdalpha - tx*sdbeta);
  if (TMath::Abs(norm) < 1E-7) { 
    WARNING("norm very small: " << norm);
  }


  // dmx / ddx
  J[0][0] = 
    -(cdbeta*cdgamma) - (sdbeta*
			 (tx*cdbeta*cdgamma + cdalpha*cdgamma*sdbeta + 
			  sdalpha*sdgamma + 
			  ty*(cdgamma*sdalpha*sdbeta - cdalpha*sdgamma)))/
    norm;
  // dmy / ddx
  J[1][0] = 
    -(cdbeta*sdgamma) - (sdbeta*
			 (-(cdgamma*sdalpha) + tx*cdbeta*sdgamma + 
			  cdalpha*sdbeta*sdgamma + 
			  ty*(cdalpha*cdgamma + sdalpha*sdbeta*sdgamma)))/
    norm;
  // dmx / ddy
  J[0][1] = 
    -(cdgamma*sdalpha*sdbeta) + cdalpha*sdgamma + 
    (cdbeta*sdalpha*(tx*cdbeta*cdgamma + 
		     cdalpha*cdgamma*sdbeta + sdalpha*sdgamma + 
		     ty*(cdgamma*sdalpha*sdbeta - cdalpha*sdgamma)))/
    norm;
  // dmy / ddy
  J[1][1] = 
    -(cdalpha*cdgamma) - sdalpha*sdbeta*sdgamma + 
    (cdbeta*sdalpha*(-(cdgamma*sdalpha) + 
		     tx*cdbeta*sdgamma + cdalpha*sdbeta*sdgamma + 
		     ty*(cdalpha*cdgamma + sdalpha*sdbeta*sdgamma)))/
    norm;
  // dmx / ddz
  J[0][2] = 
    -(cdalpha*cdgamma*sdbeta) - sdalpha*sdgamma + 
    (cdalpha*cdbeta*(tx*cdbeta*cdgamma + 
		     cdalpha*cdgamma*sdbeta + sdalpha*sdgamma + 
		     ty*(cdgamma*sdalpha*sdbeta - cdalpha*sdgamma)))/
    norm;
  // dmy / ddz
  J[1][2] = 
    cdgamma*sdalpha - cdalpha*sdbeta*sdgamma + 
    (cdalpha*cdbeta*(-(cdgamma*sdalpha) + 
		     tx*cdbeta*sdgamma + cdalpha*sdbeta*sdgamma + 
		     ty*(cdalpha*cdgamma + sdalpha*sdbeta*sdgamma)))/
    norm;
  // dmx / ddalpha
  J[0][3] = 
    -(dz*(-(cdgamma*sdalpha*sdbeta) + cdalpha*sdgamma)) + 
    (-dy + y)*(cdalpha*cdgamma*sdbeta + sdalpha*sdgamma) - 
    (((-dy + y)*cdalpha*cdbeta + dz*cdbeta*sdalpha)*
     (tx*cdbeta*cdgamma + cdalpha*cdgamma*sdbeta + 
      sdalpha*sdgamma + 
      ty*(cdgamma*sdalpha*sdbeta - cdalpha*sdgamma)))/
    norm + 
    ((ty*cdalpha*cdbeta - cdbeta*sdalpha)*
     (-(dz*cdalpha*cdbeta) + (-dy + y)*cdbeta*sdalpha - 
      (-dx + x)*sdbeta)*(tx*cdbeta*cdgamma + 
			 cdalpha*cdgamma*sdbeta + sdalpha*sdgamma + 
			 ty*(cdgamma*sdalpha*sdbeta - cdalpha*sdgamma)))/
    TMath::Power(cdalpha*cdbeta + ty*cdbeta*sdalpha - tx*sdbeta,
		 2) - ((-(dz*cdalpha*cdbeta) + 
			(-dy + y)*cdbeta*sdalpha - (-dx + x)*sdbeta)*
		       (-(cdgamma*sdalpha*sdbeta) + cdalpha*sdgamma + 
			ty*(cdalpha*cdgamma*sdbeta + sdalpha*sdgamma)))/
    norm;
  // dmy / ddalpha
  J[1][3] = 
    (-dy + y)*(-(cdgamma*sdalpha) + cdalpha*sdbeta*sdgamma) - 
    dz*(-(cdalpha*cdgamma) - sdalpha*sdbeta*sdgamma) - 
    ((-(dz*cdalpha*cdbeta) + (-dy + y)*cdbeta*sdalpha - 
      (-dx + x)*sdbeta)*(-(cdalpha*cdgamma) - 
			 sdalpha*sdbeta*sdgamma + 
			 ty*(-(cdgamma*sdalpha) + cdalpha*sdbeta*sdgamma))\
      )/norm - 
    (((-dy + y)*cdalpha*cdbeta + dz*cdbeta*sdalpha)*
     (-(cdgamma*sdalpha) + tx*cdbeta*sdgamma + 
      cdalpha*sdbeta*sdgamma + 
      ty*(cdalpha*cdgamma + sdalpha*sdbeta*sdgamma)))/
    norm + 
    ((ty*cdalpha*cdbeta - cdbeta*sdalpha)*
     (-(dz*cdalpha*cdbeta) + (-dy + y)*cdbeta*sdalpha - 
      (-dx + x)*sdbeta)*(-(cdgamma*sdalpha) + 
			 tx*cdbeta*sdgamma + cdalpha*sdbeta*sdgamma + 
			 ty*(cdalpha*cdgamma + sdalpha*sdbeta*sdgamma)))/
    TMath::Power(cdalpha*cdbeta + ty*cdbeta*sdalpha - tx*sdbeta,2);
  // dmx / ddbeta
  J[0][4] = 
    -(dz*cdalpha*cdbeta*cdgamma) + 
    (-dy + y)*cdbeta*cdgamma*sdalpha - 
    (-dx + x)*cdgamma*sdbeta - 
    ((-(dz*cdalpha*cdbeta) + (-dy + y)*cdbeta*sdalpha - 
      (-dx + x)*sdbeta)*(cdalpha*cdbeta*cdgamma + 
			 ty*cdbeta*cdgamma*sdalpha - tx*cdgamma*sdbeta))/
    norm + 
    ((-(dz*cdalpha*cdbeta) + (-dy + y)*cdbeta*sdalpha - 
      (-dx + x)*sdbeta)*(-(tx*cdbeta) - cdalpha*sdbeta - 
			 ty*sdalpha*sdbeta)*
     (tx*cdbeta*cdgamma + cdalpha*cdgamma*sdbeta + 
      sdalpha*sdgamma + 
      ty*(cdgamma*sdalpha*sdbeta - cdalpha*sdgamma)))/
    TMath::Power(cdalpha*cdbeta + ty*cdbeta*sdalpha - tx*sdbeta,
		 2) - ((-((-dx + x)*cdbeta) + dz*cdalpha*sdbeta - 
			(-dy + y)*sdalpha*sdbeta)*
		       (tx*cdbeta*cdgamma + cdalpha*cdgamma*sdbeta + 
			sdalpha*sdgamma + 
			ty*(cdgamma*sdalpha*sdbeta - cdalpha*sdgamma)))/
    norm;
  // dmy / ddbeta
  J[1][4] = 
    -(dz*cdalpha*cdbeta*sdgamma) + 
    (-dy + y)*cdbeta*sdalpha*sdgamma - 
    (-dx + x)*sdbeta*sdgamma - 
    ((-(dz*cdalpha*cdbeta) + (-dy + y)*cdbeta*sdalpha - 
      (-dx + x)*sdbeta)*(cdalpha*cdbeta*sdgamma + 
			 ty*cdbeta*sdalpha*sdgamma - tx*sdbeta*sdgamma))/
    norm + 
    ((-(dz*cdalpha*cdbeta) + (-dy + y)*cdbeta*sdalpha - 
      (-dx + x)*sdbeta)*(-(tx*cdbeta) - cdalpha*sdbeta - 
			 ty*sdalpha*sdbeta)*
     (-(cdgamma*sdalpha) + tx*cdbeta*sdgamma + 
      cdalpha*sdbeta*sdgamma + 
      ty*(cdalpha*cdgamma + sdalpha*sdbeta*sdgamma)))/
    TMath::Power(cdalpha*cdbeta + ty*cdbeta*sdalpha - tx*sdbeta,
		 2) - ((-((-dx + x)*cdbeta) + dz*cdalpha*sdbeta - 
			(-dy + y)*sdalpha*sdbeta)*
		       (-(cdgamma*sdalpha) + tx*cdbeta*sdgamma + 
			cdalpha*sdbeta*sdgamma + 
			ty*(cdalpha*cdgamma + sdalpha*sdbeta*sdgamma)))/
    norm;
  // dmx / ddgamma
  J[0][5] = 
    -((-dx + x)*cdbeta*sdgamma) - 
    dz*(cdgamma*sdalpha - cdalpha*sdbeta*sdgamma) + 
    (-dy + y)*(-(cdalpha*cdgamma) - 
	       sdalpha*sdbeta*sdgamma) - 
    ((-(dz*cdalpha*cdbeta) + (-dy + y)*cdbeta*sdalpha - 
      (-dx + x)*sdbeta)*(cdgamma*sdalpha - 
			 tx*cdbeta*sdgamma - cdalpha*sdbeta*sdgamma + 
			 ty*(-(cdalpha*cdgamma) - sdalpha*sdbeta*sdgamma)))/
    norm;
  // dmy / ddgamma
  J[1][5] =
    (-dx + x)*cdbeta*cdgamma + 
    (-dy + y)*(cdgamma*sdalpha*sdbeta - cdalpha*sdgamma) - 
    dz*(cdalpha*cdgamma*sdbeta + sdalpha*sdgamma) - 
    ((-(dz*cdalpha*cdbeta) + (-dy + y)*cdbeta*sdalpha - 
      (-dx + x)*sdbeta)*(tx*cdbeta*cdgamma + 
			 cdalpha*cdgamma*sdbeta + sdalpha*sdgamma + 
			 ty*(cdgamma*sdalpha*sdbeta - cdalpha*sdgamma)))/
    norm;

  return J;

}

/** \brief Jacobian matrix, derivatives of measurement equation to alignment
 *  parameters.
 *
 *  The matrix is evaluated at the position p0 where p0 are the track
 *  parameters, p = (tx, ty, x, y). The computation is done assuming zero
 *  misalignment as expansion point of the alignment.
 */
const HepMatrix KalmanFilterAlignmentInputProvider::Jacobian_Alignment(const HepVector & p0) const
{
  HepVector a0(6, 0);
  return Jacobian_Alignment(p0, a0);
}

