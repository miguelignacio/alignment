// ROOT includes
#include "TEnv.h"
#include "TRandom.h"

// local includes
#include "Constants.h"
#include "Detector.h"
#include "Utilities.h"
#include "ThreeDModel.h"

using namespace CLHEP;

/** \file Detector.cc
 *
 *  \brief Definition of class Detector.
 *
 */

/** Global detector handle. */
Detector * gDetector = 0;

/** Initialize Detector from given configuration file */
Detector::Detector(const char * fname)
{
  ReadConfiguration(fname);

  // reserve space for true misalignment parameters
  HepVector temp(nAlignmentParam*fNLayer*fNRod*fNModule, 0);
  fTrueMisalignment = temp;

  CreateModules();
}

/** Destructor. Delete the detectors from memory */
Detector::~Detector()
{
  for (unsigned int i = 0; i < DetVec.size(); i++) {
    delete DetVec[i];
  }
  DetVec.clear();
}

/** Read the specified configuration file, and initialize object. */
void Detector::ReadConfiguration(const char * fname)
{
  TEnv mEnv(fname);

  // Detector layout
  fNLayer           = mEnv.GetValue("NLayer", 6);
  fNRod             = mEnv.GetValue("NRod", 2);
  fNModule          = mEnv.GetValue("NModule", 6);
  fLayerDistance    = mEnv.GetValue("LayerDistance", 86);
  fRandomPlacement  = mEnv.GetValue("RandomPlacement", false);

  // Module layout
  fModuleWidth      = mEnv.GetValue("ModuleWidth", 100);
  fModuleLength     = mEnv.GetValue("ModuleLength", 200);
  fModuleThickness  = mEnv.GetValue("ModuleThickness", 0.5);
  fModulePitchX     = mEnv.GetValue("ModulePitchX", 0.183);
  fModulePitchY     = mEnv.GetValue("ModulePitchY", 1.83);

  // Misalignment
  fMisalignPosition = mEnv.GetValue("MisalignPosition", 1);
  fMisalignRotation = mEnv.GetValue("MisalignRotation", 0.050);

  INFO("Detector configuration:");
  INFO("fNLayer           : " << fNLayer);
  INFO("fNRod             : " << fNRod);
  INFO("fNModule          : " << fNModule);
  INFO("fLayerDistance    : " << fLayerDistance);
  INFO("fModuleWidth      : " << fModuleWidth);
  INFO("fModuleLength     : " << fModuleLength);
  INFO("fModuleThickness  : " << fModuleThickness);
  INFO("fModulePitchX     : " << fModulePitchX);
  INFO("fModulePitchY     : " << fModulePitchY);
  INFO("fMisalignPosition : " << fMisalignPosition);
  INFO("fMisalignRotation : " << fMisalignRotation);
}

/** Create all detector modules. Iterates over layers, rods, modules. */
void Detector::CreateModules()
{
  // reserve enough space for all detectors
  DetVec.reserve(fNLayer*fNRod*fNModule);

  // loop over all detectors
  for (int l = 0; l < fNLayer; l++) {
    for (int r = 0; r < fNRod; r++) {
      for (int m = 0; m < fNModule; m++) {
	int ID = MakeID(l, r, m);
	// Create detector at nominal position
	Det * adet = CreateDetector(ID, l, r, m);
	// misalign detector
	MisalignDetector(adet);
	// save detector in collection
	DetVec[ID] = adet;
      }
    }
  }
}

/** \brief Create one detector at the specified place.
 * 
 *  The procedure automatically ensures that the detector origin is at (0, 0,
 *  0).
 */
Det * Detector::CreateDetector(int ID, int layer, int rod, int module)
{
  Det * adet = new Det(ID, fModuleWidth, fModuleLength, fModuleThickness, fModulePitchX, fModulePitchY);
  HepVector NominalPosition;
  setXYZ(NominalPosition, 
	 rod    * fModuleWidth + fModuleWidth/2.,
	 module * fModuleLength + fModuleLength/2.,
	 layer  * fLayerDistance); //  + fModuleThickness/2.);
  ReferenceFrame nominal;
  nominal.SetPosition(NominalPosition);

  if (fRandomPlacement) {
    // DEBUG - put modules in a random place / orientation in the detector
    // This helps to check reference trajectories, track pulls etc.

    ReferenceFrame deltaFrame;
    // make a small movement
    double delta_x = gRandom->Uniform(-1, 1);
    double delta_y = gRandom->Uniform(-1, 1);
    double delta_z = gRandom->Uniform(-1, 1);
    HepVector move;
    setXYZ(move, delta_x, delta_y, delta_z);
    deltaFrame.SetPosition(move);
    // make a small rotation
    double alpha = gRandom->Uniform(-0.2, 0.2);
    double beta  = gRandom->Uniform(-0.2, 0.2);
    double gamma = gRandom->Uniform(-0.2, 0.2);
    HepMatrix rot;
    FillRotMatrixKarimaki(rot, alpha, beta, gamma);
    deltaFrame.SetRotation(rot);

    adet->SetNominalFrame(ReferenceFrame::combine_karimaki(nominal, deltaFrame));
  }
  else {
    adet->SetNominalFrame(nominal);
  }
  return adet;
}

/** Misalign the detector. Apply uniform shifts in the position and uniform
 * rotations. */
void Detector::MisalignDetector(Det * adet)
{
  ReferenceFrame deltaFrame;
  // make a small movement
  double delta_x = gRandom->Uniform(-fMisalignPosition, fMisalignPosition);
  double delta_y = gRandom->Uniform(-fMisalignPosition, fMisalignPosition);
  double delta_z = gRandom->Uniform(-fMisalignPosition, fMisalignPosition);
  HepVector move;
  setXYZ(move, delta_x, delta_y, delta_z);
  deltaFrame.SetPosition(move);
  // make a small rotation
  double alpha = gRandom->Uniform(-fMisalignRotation, fMisalignRotation);
  double beta  = gRandom->Uniform(-fMisalignRotation, fMisalignRotation);
  double gamma = gRandom->Uniform(-fMisalignRotation, fMisalignRotation);
  HepMatrix rot;
  FillRotMatrixKarimaki(rot, alpha, beta, gamma);
  deltaFrame.SetRotation(rot);
  adet->SetMisalignFrame(deltaFrame);
  // save true alignment parameters for later use
  const int offset = nAlignmentParam*adet->GetID();
  fTrueMisalignment[offset]   = delta_x;
  fTrueMisalignment[offset+1] = delta_y;
  fTrueMisalignment[offset+2] = delta_z;
  fTrueMisalignment[offset+3] = alpha;
  fTrueMisalignment[offset+4] = beta;
  fTrueMisalignment[offset+5] = gamma;
}

/** Create a single unique number used as an identity out of layer, rod and
 * module index. */
int Detector::MakeID(int layer, int rod, int module) const
{
  return layer * fNRod * fNModule + rod * fNModule + module;
}

const Det & Detector::GetDet(int layer, int rod, int module) const
{
  int ID = MakeID(layer, rod, module);
  assert(DetVec[ID] != 0);
  return *(DetVec[ID]);
}

/** Reset detector to nominal position and orientation. Additionally the true
 * alignment parameters are reset. */
void Detector::ResetDet(int moduleID)
{
  ReferenceFrame deltaFrame;
  assert(DetVec[moduleID] != 0);
  DetVec[moduleID]->SetMisalignFrame(deltaFrame);
  const int offset = nAlignmentParam*moduleID;
  fTrueMisalignment[offset]   = 0;
  fTrueMisalignment[offset+1] = 0;
  fTrueMisalignment[offset+2] = 0;
  fTrueMisalignment[offset+3] = 0;
  fTrueMisalignment[offset+4] = 0;
  fTrueMisalignment[offset+5] = 0;
}

