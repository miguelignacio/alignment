// C++ includes
#include <iostream>

// CLHEP includes
#include <CLHEP/Matrix/Vector.h>
#include <CLHEP/Matrix/Matrix.h>

// ROOT includes
#include <TCanvas.h>
#include <TViewerX3D.h>
#include <TGLViewer.h>
#include <TVirtualGL.h>
#include <TGLOutput.h>

// local includes
#include "Visualization.h"
#include "Detector.h"
#include "ThreeDModel.h"

using namespace std;
using namespace CLHEP;

/** \file Visualization.cc
 *
 * \brief Definition of class Visualization.
 *
 */

Visualization::Visualization() 
  : fGeometry("align", "align.cxx"), fHit(0), fTrack(0)
{
}

void Visualization::AddDetector(const Detector & det)
{
  double width  = det.GetWidth();
  double length = det.GetLength();
  double height = det.GetHeight();

  /// Create Materials
  TMaterial * mat;
  mat = new TMaterial("silicon","SILICON",28.09,14,2.329999);
  mat = new TMaterial("air","AIR",14.60999,7.3,.001205);

  /// Create world volume
  TBRIK * brik;
  brik = new TBRIK("SHEBANG", "SHEBANG", "air", 
		   2*width, 2*length, 2*height);
  brik->SetVisibility(0);

  /// enter world
  TNode * node = new TNode("WORLD", "WORLD", "SHEBANG");

  // create placeholder for a hit
  fHit = new TBRIK("HIT", "HIT", "air", 5, 5, 5);

  node->cd();

  // now add all detectors
  int nlayer  = det.GetNLayer();
  int nrod    = det.GetNRod();
  int nmodule = det.GetNModule();
  for (int l = 0; l < nlayer; l++) {
    for (int r = 0; r < nrod; r++) {
      for (int m = 0; m < nmodule; m++) {
	AddDet(det.GetDet(l, r, m), l+2);
      }
    }
  }
}

void Visualization::AddDet(const Det & det, int color)
{
  double width  = det.GetWidth();
  double length = det.GetLength();
  double depth  = det.GetThickness();

  // create placeholders for a detector
  TBRIK * brik = new TBRIK("DET", "DET", "silicon", width/2., length/2., depth/2.);

  const ReferenceFrame & frame = det.GetMisalign();
  HepMatrix Rotation = frame.GetRotation();
  double matrix_values[9];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      matrix_values[i*3+j] = Rotation[i][j];
    }
  }

  TRotMatrix * rot;
  rot = new TRotMatrix("rot", "rot",
// 		       Form("rot%d_%d_%d", l, r, m),
// 		       Form("rot%d_%d_%d", l, r, m),
		       matrix_values);

  // total translation = misalignment translation + nominal translation
  TNode * node;
  HepVector position = frame.GetPosition();
  node = new TNode("det", "det", 
// 		   Form("Det%02d_%02d_%02d", l, r, m),
// 		   Form("Det%02d_%02d_%02d", l, r, m),
		   brik, 
		   position[0], position[1], position[2], 
		   rot);
// 		   Form("rot%d_%d_%d", l, r, m));
  node->SetLineColor(color);
}

void Visualization::AddHit(const Hit & hit)
{
  const Det & det = hit.GetDet();
  HepVector true_pos = det.GetMisalign().ToGlobal(hit.GetPosition());
  HepVector measured_pos = det.GetNominal().ToGlobal(hit.GetPosition());
  TNode * n;
  n = new TNode("h", "h", fHit, true_pos[0], true_pos[1], true_pos[2]);
//   cout << "Hit in layer " << hit.GetDet()->layer() << endl;
  n->SetLineColor(3);
  n = new TNode("h", "h", fHit, measured_pos[0], measured_pos[1], measured_pos[2]);
  n->SetLineColor(1);
}

void Visualization::AddTrack(const Track & track, const int color)
{
  /// store number of track to have unique names for TRotMatrix...
//   static int track_number = 0;
//   track_number++;

  // track direction
  HepVector direction = track.GetDestination() - track.GetOrigin();
//   cout << "direction vector is " << direction << endl;

  // track length
  double length = sqrt(direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2]);
//   cout << "length should be " << length << endl;

  // calculate polar angle theta and azimuth angle phi of track
  double theta = acos(direction[2]/length);
  double phi   = atan2(direction[1], direction[0]);
//   cout << "track theta = " << theta << endl;
//   cout << "track phi = " << phi << endl;

  // fill a rotation matrix with polar angle theta and azimuth angle phi
  HepMatrix m;
  FillRotMatrixPolar(m, theta, phi);
  // convert it such that TRotMatrix can handle it
  double matrix_values[9];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      matrix_values[i*3+j] = m[i][j];
    }
  }
//   cout << "rotation matrix " << m << endl;

  // fill TRotMatrix (needed for display in Geometry)
  TRotMatrix * rot;
  rot = new TRotMatrix("rot", "rot",
//                     Form("rot%d", track_number), 
// 		       Form("rot%d", track_number), 
		       matrix_values);
//   cout << "origin " << origin << "destination " << destination << endl;
  HepVector middle = track.GetOrigin() + 0.5*direction;
//   cout << "middle is here " << middle << endl;

  // create track
  fTrack = new TBRIK("TRACK", "TRACK", "air", 1., 1., length/1.9);
  // create node with position and rotation, showing the track
  TNode * n;
  n = new TNode("tr", "tr", fTrack, middle[0], middle[1], middle[2], rot);
// 		Form("rot%d", track_number));
  n->SetLineColor(color);
}

void Visualization::Show()
{
  // create viewer
  TCanvas * c1;
  c1 = new TCanvas("c1", "c1", 0, 0, 800, 800);
  fGeometry.Draw();
  TGLViewer * oglv = (TGLViewer *) gPad->GetViewer3D("ogl");
  TGLViewer::ECameraType ct = TGLViewer::kCameraPerspXOY;
  oglv->SetCurrentCamera(ct);
  double fov = 40;
  double dolly = 200;
  double center[3] = {gDetector->GetWidth()/2, gDetector->GetLength()/2., gDetector->GetHeight()/2.};
  double hRotate = -0.2;
  double vRotate = -0.5;
  oglv->SetPerspectiveCamera(ct, fov, dolly, center, hRotate, vRotate);
  // TViewerX3D * viewer = new TViewerX3D(gPad);
  // viewer->DrawViewer();
}
