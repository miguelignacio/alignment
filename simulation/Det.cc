#include <iostream>

#include "Det.h"
#include "Utilities.h"
#include "ThreeDModel.h"

using namespace std;
using namespace CLHEP;

/** \file Det.cc  \brief Definition of class Det. */

/** \brief Recommended constructor.
 *
 * You should use this constructor only.
 */
Det::Det(int ID, double Width, double Length, double Thickness, double PitchX, double PitchY)
  : fID(ID), fWidth(Width), fLength(Length), fThickness(Thickness), fPitchX(PitchX), fPitchY(PitchY)
{ 
}

/** \brief Get radiation length (in cm)
 *
 *  Material is assumed to be silicon.
 */
const double Det::GetRadiationLength() const
{
  const double si_xi          = 21.82; // g/cm^2
  const double si_rho         = 2.33; // g/cm^3
  return si_xi/si_rho;
}

/** \brief Set nominal frame.
 *
 * This also recalculates the misaligned and aligned frame corresponding to
 * the given delta frames.
 *
 */
void Det::SetNominalFrame(ReferenceFrame nominal) 
{ 
  fNominal  = nominal; 
  fMisalign = ReferenceFrame::combine_karimaki(fNominal, fMisalignDelta);
  fAlign    = ReferenceFrame::combine_karimaki(fNominal, fAlignDelta); 
}

/** \brief Set misaligned frame.
 *
 * This sets both the delta (that you supply) and combines it with the nominal
 * frame.
 *
 */
void Det::SetMisalignFrame(ReferenceFrame misalignDelta)
{ 
  fMisalignDelta = misalignDelta;
  fMisalign = ReferenceFrame::combine_karimaki(fNominal, fMisalignDelta);
}

/** \brief Set aligned frame.
 *
 * This sets both the delta (that you supply) and combines it with the nominal
 * frame.
 *
 */
void Det::SetAlignFrame(ReferenceFrame alignDelta)
{ 
  fAlignDelta = alignDelta; 
  fAlign = ReferenceFrame::combine_karimaki(fNominal, fAlignDelta); 
}

/** \brief Check if the given hit (in local coordinates) is inside detector
 * volume.
 *
 */
bool Det::InsideBounds(HepVector local) const
{
  //   cout << "local position is " << local << endl;
  if (local[0] < -fWidth/2.  || local[0] > fWidth/2.) {
    // cout << "local[0] = " << local[0] << ", fWidth/2. = " << fWidth/2. << endl;
    return false;
  }
  if (local[1] < -fLength/2. || local[1] > fLength/2.) {
    // cout << "local[1] = " << local[1] << ", fLength/2. = " << fLength/2. << endl;
    return false;
  }
  if (local[2] < -fThickness/2. || local[2] > fThickness/2.) {
    // cout << "local[2] = " << local[2] << ", fThickness/2. = " << fThickness/2. << endl;
    return false;
  }
  return true;
}

/** \brief Get measured intersection assuming nominal position for detector.
 *
 * Get from real misaligned state the intersection that would be measured if
 * the misaligned state is not known, i.e. regarding only nominal position /
 * rotation.
 *
 */
HepVector Det::MeasuredIntersection(HepVector TrueIntersection) const
{
  // get true local position
  HepVector local = GetMisalign().ToLocal(TrueIntersection);
//   cout << "Local coordinates " << local << endl;
  // convert (wrongly) to global coordinates
  HepVector global = GetNominal().ToGlobal(local);
//   cout << "true intersection was " << TrueIntersection << endl;
//   cout << "Difference is " << global - TrueIntersection;
  return global;
}
