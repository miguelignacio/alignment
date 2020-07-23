#ifndef _Track_h
#define _Track_h

#include <vector>

#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Matrix/SymMatrix.h>

#include "Det.h"
#include "Hit.h"

/** \file Track.h
 *
 * \brief Declaration of class Track
 *
 */

/** \class Track
 *
 * \brief Description of a straight track.
 *
 *
 * A track can be described with several parameters. Here, a start and end
 * point parametrization is used in the constructor, and internally a
 * parametrization by the definition of a (x,y) point on some reference plane
 * and the tangent to the track is used as parameters.
 *
 */

class Track {
protected:
  // track parameters on a given surface
  /// track origin
  CLHEP::HepVector    fOrigin;
  /// track destination
  CLHEP::HepVector    fDestination;
  /// track covariance matrix in global frame
  CLHEP::HepSymMatrix fCovariance;
  /// Track momentum in GeV. It is used only for multiple scattering. The
  /// particle is assumed to be a muon.
  double fMomentum;
  /// track chi^2. Initialized to -1, only meaningful once track reconstructed
  double fChi2;
  /// number of degrees of freedom. Initialized to zero.
  int fNdof;
  /// vector of hits that the track has in the detector. Can be empty and
  /// needs to be filled externally, e.g. during simulation or
  /// reconstruction. However, the track takes ownership of the hits,
  /// i.e. deletes them in the destructor.
  std::vector<Hit *> fHits;

public:
  Track();
  Track(CLHEP::HepVector origin, CLHEP::HepVector destination, double momentum);
  virtual ~Track();
  /// Compute intersection with plane fixed at (0,0,0) and normal (0,0,1)
  CLHEP::HepVector GetOrigin() const;
  /// Compute intersection with plane fixed at (0,0,detector height) and normal (0,0,1).
  CLHEP::HepVector GetDestination() const;
  /// Get particle momentum (in GeV)
  double GetMomentum() const { return fMomentum; }
  /// Get particle mass (in GeV)
  double GetMass() const { return 0.105685; }
  /// Get particle charge
  int GetCharge() const { return 1; }
  /// Get chi2
  double GetChi2() const { return fChi2; }
  /// Get ndof
  int GetNdof() const { return fNdof; }
  /// Set chi2 and ndof
  void SetChi2Ndof(double chi2, int ndof) { fChi2 = chi2; fNdof = ndof; }
  /// Compute the intersection with the detector and return a Hit. Boundary
  /// checking is performed. In case the hit is outside the detector bounds, a
  /// null pointer is retuned.
  Hit * MisalignedIntersection(const Det & adet) const;
  Hit * NominalIntersection(const Det & adet) const;
  /// Get intersection of track with a plane in global coordinate system.
  CLHEP::HepVector * IntersectionWithPlane(const ReferenceFrame & frame) const;
  // Parameters
  // format is (tx, ty, x, y) 
  // where tx = dx/dz = tan(theta)*cos(phi)
  // where ty = dy/dz = tan(theta)*sin(phi)
  // x, y global x and y coordinates at z=0
  CLHEP::HepVector GetGlobalParameters() const;
  void SetGlobalParameters(CLHEP::HepVector p);
  void SetGlobalCovariance(CLHEP::HepSymMatrix covariance) { fCovariance = covariance; };
  CLHEP::HepSymMatrix GetGlobalCovariance() const { return fCovariance; }
  CLHEP::HepVector GetLocalParameters(const ReferenceFrame & frame) const;
  CLHEP::HepSymMatrix GetLocalCovariance(const ReferenceFrame & frame);
  void SetLocalCovariance(const ReferenceFrame & frame, const CLHEP::HepSymMatrix & covariance);
  // derivatives from local to global parameters
  CLHEP::HepMatrix Jacobian(const ReferenceFrame & frame) const;
  void SetLocalParameters(const ReferenceFrame & frame, CLHEP::HepVector p);
  const std::vector<Hit *> & GetHits() const { return fHits; };
  std::vector<Hit *> & GetHits() { return fHits; };
};

#endif
