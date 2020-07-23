#include <CLHEP/Matrix/DiagMatrix.h>

#include "ReferenceFrame.h"

using namespace CLHEP;

/** \file ReferenceFrame.h
 *
 * \brief Definition of class ReferenceFrame
 *
 */

// Default constructor
ReferenceFrame::ReferenceFrame() : fPosition(3, 0), fRotation(3, 3, 1)
{
}

ReferenceFrame::ReferenceFrame(HepVector position, HepMatrix rotation) : fPosition(position), fRotation(rotation)
{
}


ReferenceFrame ReferenceFrame::combine_karimaki(const ReferenceFrame & first, const ReferenceFrame & delta)
{
  ReferenceFrame combinedFrame;
  combinedFrame.SetRotation(delta.GetRotation()*first.GetRotation());
  combinedFrame.SetPosition(delta.GetPosition()+first.GetPosition());
  return combinedFrame;
}
