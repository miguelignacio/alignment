#ifndef _Det_h
#define _Det_h

#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>

#include "ReferenceFrame.h"

/** \file Det.h
 *
 *  \brief Declaration of class Det
 *
 */

/** \class Det
 *
 * \brief A (misaligned) detector in space. 
 *
 *  It consists of three locations: The nominal position of the object (where
 *  it is expected to be), the true (generated) misalignment, and the computed
 *  alignment. The two members fMisalignDelta and fAlignDelta are no real
 *  reference frames, but rather store additional rotation and position that
 *  have to be applied additionally to the nominal to get the misaligned and
 *  aligned frame, respectively.
 *
 *  The two members fMisalign and fAlign contain the combination done in the
 *  Karimaki way of the nominal plus the applied delta. They are updated each
 *  time SetMisalignDeltaFrame or SetAlignDeltaFrame is called.
 */
class Det {
 protected:
  /// internal number of module, module "ID"
  int fID; 
  /// Detector width (x-coordinate)
  double fWidth; 
  /// Detector length (y-coordinate)
  double fLength;
  /// Detector thickness (z-coordinate)
  double fThickness;
  /// Pitch in x coordinate (along "width")
  double fPitchX;
  /// Pitch in y coordinate (along "length")
  double fPitchY;
  ReferenceFrame fNominal;
  ReferenceFrame fMisalignDelta;
  ReferenceFrame fAlignDelta;
  ReferenceFrame fMisalign;
  ReferenceFrame fAlign;

 public:
  Det(int ID, double Width, double Length, double Thickness, double PitchX, double PitchY);

  // return module parameters

  /// width of a module (perpendicular to strips, local x)
  const double GetWidth() const { return fWidth; };
  /// length of a module (along strips, local y)
  const double GetLength() const { return fLength; };
  /// thickness of a module (local z)
  const double GetThickness() const { return fThickness; };
  /// the strip pitch of a module (sensitive coordinate along x)
  const double GetPitchX() const { return fPitchX; };
  /// the stereo coordinate distance, i.e. y distance of crossings
  const double GetPitchY() const { return fPitchY; };
  /// detector resolution along x
  const double GetResolutionX() const { return fPitchX/sqrt(12.); };
  /// detector resolution along y
  const double GetResolutionY() const { return fPitchY/sqrt(12.); };
  /// detector radiation length in cm
  const double GetRadiationLength() const;
  /// Detector ID
  int GetID() const { return fID; };

  // Access the reference frames. Note that here not only the delta frames can
  // be returned, but also the combined frames.
  /// Return the nominal frame (i.e. where the detector is supposed to be)
  const ReferenceFrame & GetNominal() const { return fNominal; };
  /// Return difference from nominal frame to real (misaligned frame)
  const ReferenceFrame & GetMisalignDelta() const { return fMisalignDelta; };
  /// Return absolute misaligned frame (nominal position + delta due to misalignment)
  const ReferenceFrame & GetMisalign() const { return fMisalign; };
  /// Return difference from nominal to aligned frame
  const ReferenceFrame & GetAlignDelta() const { return fAlignDelta; };
  /// Return absolute aligned frame (nominal position + computed alignment correction)
  const ReferenceFrame & GetAlign() const { return fAlign; };

  /// Boundary checking: Is a point in the detector?
  /// returns true if the point is in the detector, else false
  bool InsideBounds(CLHEP::HepVector local) const;

  /// Convert the true intersection in the one that would be measured
  /// because the misalignment is not known.
  CLHEP::HepVector MeasuredIntersection(CLHEP::HepVector trueIntersection) const;

  /// Set nominal frame (where the detector should be, ideal geometry)
  void SetNominalFrame(ReferenceFrame nominal);
  /// Set misaligned frame (You have to give the DELTA to the nominal frame)
  void SetMisalignFrame(ReferenceFrame misalignDelta);
  /// Set aligned frame (You have to give the DELTA to the nominal frame)
  void SetAlignFrame(ReferenceFrame alignDelta);
};

#endif
