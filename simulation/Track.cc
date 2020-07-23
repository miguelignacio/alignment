#include <iostream>

#include <CLHEP/Matrix/Vector.h>

#include "Detector.h"
#include "Utilities.h"
#include "ThreeDModel.h"
#include "Track.h"
#include "TMath.h"

using namespace std;
using namespace CLHEP;

/** \file Track.cc
 *
 * \brief Decfinition of class Track
 *
 */

/** A track parallel to z in the center of the detector, with a momentum
 *  according to the mean of cosmic muons at the surface of the earth.
 */
Track::Track() : fCovariance(4, 0)
{
  assert(gDetector != 0);
  double width  = gDetector->GetWidth();
  double length = gDetector->GetLength();
  double height = gDetector->GetHeight();
  setXYZ(fOrigin, width/2., length/2., 0);
  HepVector destination;
  setXYZ(fDestination, width/2., length/2., height);
  fMomentum = 4;
  fChi2 = -1;
  fNdof = 0;
}

/** A track specified by the user, going from origin to destination, with the
 *  given momentum. */
Track::Track(HepVector origin, HepVector destination, double momentum)
{
  fOrigin = origin;
  fDestination = destination;
  fMomentum = momentum;
  fChi2 = -1;
  fNdof = 0;
}

/** Destruct the track. Has to delete all the hits. */
Track::~Track()
{
  for (vector<Hit *>::const_iterator it = fHits.begin(); it != fHits.end(); it++) {
    delete (*it);
  }
}

/// Compute intersection with plane fixed at (0,0,0) and normal (0,0,1)
HepVector Track::GetOrigin() const
{
  // origin at (0, 0, 0)
  HepVector position;
  setXYZ(position, 0, 0, 0);

  // local frame = global frame
  HepMatrix rotation(3, 3, 1); 
  ReferenceFrame frame(position, rotation);

  HepVector * intersection = IntersectionWithPlane(frame);
  if (intersection == 0) 
    THROW("Track parallel to global x-y plane -- no origin!");
  
  HepVector tmp = *intersection;
  delete intersection;
  return tmp;
}

/// Compute intersection with plane fixed at (0,0,detector height) and normal (0,0,1).
HepVector Track::GetDestination() const
{
  assert(gDetector != 0);
  double height = gDetector->GetHeight();

  // origin at (0, 0, height)
  HepVector position;
  setXYZ(position, 0, 0, height);

  // local frame = global frame
  HepMatrix rotation(3, 3, 1); 
  ReferenceFrame frame(position, rotation);

  HepVector * intersection = IntersectionWithPlane(frame);
  if (intersection == 0)
    THROW("Track parallel to global x-y plane -- no destination!");
  
  HepVector tmp = *intersection;
  delete intersection;
  return tmp;
}

/** \brief Compute the intersection point of the track with plane given by
 * "reference frame" in global coordinates.
 *
 *  The plane is parametrized by \vec{x} * \vec{n} = d.
 *  The track is parametrized by \vec{x} = \vec{x}_0 + \lambda * (\vec{x}_1 - \vec{x}_0) .
 *  No boundary checking is performed.
 */
HepVector * Track::IntersectionWithPlane(const ReferenceFrame & frame) const
{
  LOG(5, "Frame origin is " << frame.GetPosition());
  LOG(5, "Frame rotation is " << frame.GetRotation());
  LOG(5, "Track origin is" << fOrigin);
  LOG(5, "Track destination is" << fDestination);
  HepVector z;
  setXYZ(z, 0, 0, 1);
  HepVector n = frame.GetRotation().T() * z;
  /// direction of track
  HepVector m =  fDestination - fOrigin;
  LOG(5, "Normal of detector is " << n);
  /// get distance from origin of detector
  double d = dot(frame.GetPosition(), n);
  LOG(5, "Distance of plane d = " << d);
  double denominator = dot(m,n);
  if (denominator == 0) {
    // track parallel to plane, no solution
    return 0;
  }
  double lambda = (d - dot(fOrigin,n)) / denominator;
  LOG(5, "lambda = " << lambda);
  /// intersection
  HepVector * i = new HepVector(fOrigin + (lambda * m));
  LOG(5, "Intersection point is " << *i);
  return i;
}

/** \brief Compute the hit position of the track with a detector in misaligned
 * state.
 *  
 * Returns null pointer if no hit in this detector, else returns the hit.
 *
 */
Hit * Track::MisalignedIntersection(const Det & adet) const
{
  /// Compute intersection with detector plane
  HepVector * i = IntersectionWithPlane(adet.GetMisalign());
  if (i == 0)
    return 0;
  /// check detector bounds
  HepVector local = adet.GetMisalign().ToLocal(*i);
  delete i;
  if (adet.InsideBounds(local)) {
    return new Hit(adet, local);
  }
  else {
    return 0;
  }
}

/** \brief Compute the hit position of the track with a detector in nominal
 * state.
 *  
 * Returns null pointer if no hit in this detector, else returns the hit.
 *
 */
Hit * Track::NominalIntersection(const Det & adet) const
{
  /// Compute intersection with detector plane
  HepVector * i = IntersectionWithPlane(adet.GetNominal());
  if (i == 0)
    return 0;

  /// check detector bounds
  HepVector local = adet.GetNominal().ToLocal(*i);
  delete i;
  if (adet.InsideBounds(local)) {
//     cout << "inside det" << endl;
    return new Hit(adet, local);
  }
  else
    return 0;
}

/** Starting from origin and destination, fill track parameters. For this, the
 *  reference surface is z=0. 
 */
HepVector Track::GetGlobalParameters() const
{
  // get direction vector
  HepVector direction = fDestination - fOrigin;
  // new track for four parameters
  HepVector t(4);
  if (direction[2] == 0.) {
    // \warning
    // we have a problem here since the tangent then is ill-defined.
    // we cannot accept such tracks
    // by construction anyway this should not happen
    THROW("ERR: Internal misconception in Track::GetGlobalParameters()");
  }
  // tx
  t[0] = direction[0]/direction[2];
  // ty
  t[1] = direction[1]/direction[2];
  // compute intersection with z=0 plane
  double lambda = -fOrigin[2]/direction[2];
  // x
  t[2] = fOrigin[0] + lambda * direction[0];
  // y
  t[3] = fOrigin[1] + lambda * direction[1];
  return t;
}

/** Given parameters p, compute origin, direction and finally destination.
 */
void Track::SetGlobalParameters(HepVector p)
{
  assert(gDetector != 0);
  double height = gDetector->GetHeight();
  setXYZ(fOrigin, p[2], p[3], 0.);
  HepVector direction(3);
  setXYZ(direction, p[0], p[1], 1.);
  // reference plane on last detector
  double lambda = height;
  fDestination = fOrigin + lambda * direction;
}

/** Get local track parameters for given detector, using misaligned state. No
 *  boundary checking.
 */
HepVector Track::GetLocalParameters(const ReferenceFrame & frame) const
{
  /// for storing the results
  HepVector p(4);
  /// first part: get local intersection coordinates
  HepVector * hit = IntersectionWithPlane(frame);
  if (hit == 0) {
    // plane parallel to track
    THROW("ERR: Internal misconception in Track::GetLocalParameters()");
  }
  HepVector local = frame.ToLocal(*hit);
  delete hit;
  p[2] = local[0];
  p[3] = local[1];
  HepVector direction = fDestination - fOrigin;
  HepVector n = frame.GetRotation() * direction;
  if (n[2] == 0.) {
    /// track is parallel to plane, should not happen here because
    /// checked previously
    THROW("ERR: Internal error in Track::GetLocalParameters()");
  }
  p[0] = n[0] / n[2];
  p[1] = n[1] / n[2];
  return p;
}

void Track::SetLocalParameters(const ReferenceFrame & frame, HepVector p)
{
  // construct local direction vector
  HepVector n(3);
  setXYZ(n, p[0], p[1], 1.);
  // turn vector into global coordinate System
  HepMatrix Rotation = frame.GetRotation();
  HepVector direction = Rotation.T() * n;
//   cout << "direction vector is " << direction << endl;
  // build position
  HepVector pos(3);
  setXYZ(pos, p[2], p[3], 0.);
//   cout << "local position vector is " << pos << endl;
  // transform to global
  HepVector g = frame.ToGlobal(pos);
//   cout << "global position vector is " << g << endl;
  // intercept track at reference plane z=0
  if (direction[2] == 0.) {
    THROW("ERR: Track parallel to z=0 plane in Track::SetLocalParameters()");
  }
  double lambda = -g[2] / direction[2];
//   cout << "lambda is " << lambda << endl;
  setXYZ(fOrigin, g[0]+lambda*direction[0], g[1]+lambda*direction[1], 0.);

  // get config
  assert(gDetector != 0);
  double height = gDetector->GetHeight();

  // extrapolate track to last layer
  lambda = height / direction[2];
  fDestination = fOrigin + lambda * direction;
}

HepSymMatrix Track::GetLocalCovariance(const ReferenceFrame & frame)
{
  HepMatrix J = Jacobian(frame);
  return fCovariance.similarity(J);
}

void Track::SetLocalCovariance(const ReferenceFrame & frame, const HepSymMatrix & covariance)
{
  HepMatrix J = Jacobian(frame);
  int ierr;
  HepMatrix iJ = J.inverse(ierr);
  if (ierr != 0) {
    THROW("matrix inversion failed");
  }
  fCovariance = covariance.similarity(iJ);
}

/** \brief Jacobian matrix of local to global track parameters. 
 *
 *   This method calculates the derivatives of local track parameters to
 *   global track parameters for the given reference frame.
 */
HepMatrix Track::Jacobian(const ReferenceFrame & frame) const
{
  // prepare track parameters
  HepVector p0 = GetGlobalParameters();
  const double tx = p0[0];
  const double ty = p0[1];
  const double x  = p0[2];
  const double y  = p0[3];

  // prepare reference frame parameters

  // position
  HepVector position = frame.GetPosition();
  const double x0 = position[0];
  const double y0 = position[1];
  const double z0 = position[2];

  // rotation
  double alpha, beta, gamma;
  GetAnglesKarimaki(frame.GetRotation(), alpha, beta, gamma);
  const double salpha = TMath::Sin(alpha);
  const double calpha = TMath::Cos(alpha);
  const double sbeta  = TMath::Sin(beta);
  const double cbeta  = TMath::Cos(beta);
  const double sgamma = TMath::Sin(gamma);
  const double cgamma = TMath::Cos(gamma);

  // some abbreviations to save computation time
  double norm = calpha*cbeta + ty*cbeta*salpha - tx*sbeta;
  double norm2 = norm*norm;
  double a = cgamma*salpha*sbeta - calpha*sgamma;

  HepMatrix j(4,4);

  j[0][0] = 
    (cbeta*cgamma)/
    (norm) + 
    (sbeta*(tx*cbeta*cgamma + calpha*cgamma*sbeta + 
	    salpha*sgamma + 
	    ty*(a)))/
    norm2;

  j[0][1] = 
    (a)/
    (norm) - 
    (cbeta*salpha*(tx*cbeta*cgamma + 
		   calpha*cgamma*sbeta + salpha*sgamma + 
		   ty*(a)))/
    norm2;

  j[0][2] = 
    0;

  j[0][3] = 
    0;

  j[1][0] = 
    (cbeta*sgamma)/
    (norm) + 
    (sbeta*(-(cgamma*salpha) + tx*cbeta*sgamma + 
	    calpha*sbeta*sgamma + 
	    ty*(calpha*cgamma + salpha*sbeta*sgamma)))/
    norm2;

  j[1][1] = 
    (calpha*cgamma + salpha*sbeta*sgamma)/
    (norm) - 
    (cbeta*salpha*(-(cgamma*salpha) + 
		   tx*cbeta*sgamma + calpha*sbeta*sgamma + 
		   ty*(calpha*cgamma + salpha*sbeta*sgamma)))/
    norm2;

  j[1][2] = 
    0;

  j[1][3] = 
    0;

  j[2][0] = 
    cbeta*cgamma*((tx*sbeta*
		   (z0*calpha*cbeta + (-y + y0)*cbeta*salpha - 
		    (-x + x0)*sbeta))/
		  norm2 \
		  + (z0*calpha*cbeta + (-y + y0)*cbeta*salpha - 
		     (-x + x0)*sbeta)/
		  (norm)) + 
    (ty*sbeta*(z0*calpha*cbeta + 
	       (-y + y0)*cbeta*salpha - (-x + x0)*sbeta)*
     (a))/
    norm2 + 
    (sbeta*(z0*calpha*cbeta + (-y + y0)*cbeta*salpha - 
	    (-x + x0)*sbeta)*(calpha*cgamma*sbeta + 
			      salpha*sgamma))/
    norm2;

  j[2][1] = 
    -((tx*TMath::Power(cbeta,2)*cgamma*salpha*
       (z0*calpha*cbeta + (-y + y0)*cbeta*salpha - 
	(-x + x0)*sbeta))/
      norm2) \
    + (-((ty*cbeta*salpha*
          (z0*calpha*cbeta + (-y + y0)*cbeta*salpha - 
	   (-x + x0)*sbeta))/
	 TMath::Power(calpha*cbeta + ty*cbeta*salpha - 
		      tx*sbeta,2)) + (z0*calpha*cbeta + 
				      (-y + y0)*cbeta*salpha - (-x + x0)*sbeta)/
       (norm))*
    (a) - 
    (cbeta*salpha*(z0*calpha*cbeta + 
		   (-y + y0)*cbeta*salpha - (-x + x0)*sbeta)*
     (calpha*cgamma*sbeta + salpha*sgamma))/
    norm2;

  j[2][2] = 
    cbeta*cgamma*(1 + (tx*sbeta)/
		  (norm)) + 
    (ty*sbeta*(a))/
    (norm) + 
    (sbeta*(calpha*cgamma*sbeta + salpha*sgamma))/
    (norm);

  j[2][3] = 
    -((tx*TMath::Power(cbeta,2)*cgamma*salpha)/
      (norm)) + 
    (1 - (ty*cbeta*salpha)/
     (norm))*
    (a) - 
    (cbeta*salpha*(calpha*cgamma*sbeta + 
		   salpha*sgamma))/
    (norm);

  j[3][0] = 
    cbeta*((tx*sbeta*(z0*calpha*cbeta + 
		      (-y + y0)*cbeta*salpha - (-x + x0)*sbeta))/
	   TMath::Power(norm,
			2) + (z0*calpha*cbeta + (-y + y0)*cbeta*salpha - 
			      (-x + x0)*sbeta)/
	   (norm))*
    sgamma + (sbeta*(z0*calpha*cbeta + 
		     (-y + y0)*cbeta*salpha - (-x + x0)*sbeta)*
	      (-(cgamma*salpha) + calpha*sbeta*sgamma))/
    norm2 + 
    (ty*sbeta*(z0*calpha*cbeta + (-y + y0)*cbeta*salpha - 
	       (-x + x0)*sbeta)*(calpha*cgamma + 
				 salpha*sbeta*sgamma))/
    norm2;

  j[3][1] = 
    -((tx*TMath::Power(cbeta,2)*salpha*
       (z0*calpha*cbeta + (-y + y0)*cbeta*salpha - 
	(-x + x0)*sbeta)*sgamma)/
      norm2) \
    - (cbeta*salpha*(z0*calpha*cbeta + 
		     (-y + y0)*cbeta*salpha - (-x + x0)*sbeta)*
       (-(cgamma*salpha) + calpha*sbeta*sgamma))/
    norm2 + 
    (-((ty*cbeta*salpha*
	(z0*calpha*cbeta + (-y + y0)*cbeta*salpha - 
	 (-x + x0)*sbeta))/
       TMath::Power(norm,
		    2)) + (z0*calpha*cbeta + (-y + y0)*cbeta*salpha - 
			   (-x + x0)*sbeta)/
     (norm))*
    (calpha*cgamma + salpha*sbeta*sgamma);

  j[3][2] = 
    cbeta*(1 + (tx*sbeta)/
	   (norm))*
    sgamma + (sbeta*(-(cgamma*salpha) + 
		     calpha*sbeta*sgamma))/
    (norm) + 
    (ty*sbeta*(calpha*cgamma + salpha*sbeta*sgamma))/
    (norm);

  j[3][3] = 
    -((tx*TMath::Power(cbeta,2)*salpha*sgamma)/
      (norm)) - 
    (cbeta*salpha*(-(cgamma*salpha) + 
		   calpha*sbeta*sgamma))/
    (norm) + 
    (1 - (ty*cbeta*salpha)/
     (norm))*
    (calpha*cgamma + salpha*sbeta*sgamma);

  LOG(6, "Jacobian for track parameter conversion" << j);
  return j;
}
