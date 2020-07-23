#ifndef _ReferenceFrame_h
#define _ReferenceFrame_h

#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>

/** \file ReferenceFrame.h
 *
 * \brief Declaration of class ReferenceFrame
 *
 */

/** \class ReferenceFrame
 *  \brief Enables the use of local reference frames. 
 *
 * It consists of a vector and a rotation matrix. The vector is pointing to
 * the position in space of the reference frame, and the rotation is the
 * conversion from the global frame into the local frame. The class offers
 * methods to convert a position in space from the global to this reference
 * frame. No assumption about Euler angles is made since the full rotation
 * matrix is stored.
 *
 */

class ReferenceFrame {
  CLHEP::HepVector fPosition;
  CLHEP::HepMatrix fRotation;

 public:
  // constructor & destructor
  ReferenceFrame();
  ReferenceFrame(CLHEP::HepVector position, CLHEP::HepMatrix rotation);

  // methods for access to member variables
  CLHEP::HepVector GetPosition() const { return fPosition; };
  CLHEP::HepMatrix GetRotation() const { return fRotation; };
  void SetPosition(CLHEP::HepVector position) { fPosition = position; };
  void SetRotation(CLHEP::HepMatrix rotation) { fRotation = rotation; };
  CLHEP::HepVector & GetPosRef() { return fPosition; };
  CLHEP::HepMatrix & GetRotRef() { return fRotation; };

  // coordinate transformation from global to reference frame
  CLHEP::HepVector ToLocal(CLHEP::HepVector global) const { return fRotation * (global - fPosition); };
  // coordinate transformation from reference to global frame
  CLHEP::HepVector ToGlobal(CLHEP::HepVector local) const { return fRotation.T() * local + fPosition; };

  // static member to combine two reference frames. Note that the order is
  // important. Note that also it does *NOT* correspond to two consecutive
  // transformations from global to first frame, and then to the delta
  // frame. Instead, the two Rotations are multiplied and the two positions
  // added.  This is done how Veikko Karimaki defined the "misalignment". It
  // has the advantage that once you are satisfied with the combined values,
  // you can make them to the nominal ones...
  static ReferenceFrame combine_karimaki(const ReferenceFrame & first, const ReferenceFrame & delta);
};

#endif
