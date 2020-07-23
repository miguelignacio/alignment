#ifndef _Visualization_h
#define _Visualization_h

// geometry
#include <TGeometry.h>
#include <TNode.h>
#include <TMaterial.h>
#include <TBRIK.h>
#include <TTUBE.h>


// local includes
#include "ReferenceFrame.h"
#include "Det.h"
#include "Hit.h"
#include "Track.h"
#include "Detector.h"

/** \file Visualization.h
 *
 * \brief Declaration of class Visualization.
 *
 */

/** \class Visualization
 *
 *  \brief Graphical display of detector, hits and tracks.
 *
 *  Allows one to visualize the detector, and optionally show hits and tracks
 *  can be added to the graphical representation. Both x3d and OpenGL do work.
 *
 */

class Visualization
{
  TGeometry fGeometry;
  TBRIK * fHit;
  TBRIK * fTrack;

 public:
  // standard constructor
  Visualization();
  void AddDetector(const Detector & det);
  void AddHit(const Hit & hit);
  void AddTrack(const Track & track, const int color = 1);
  void Show();
 protected:
  void AddDet(const Det & det, int color);
};

#endif
