#ifndef _Detector_h
#define _Detector_h

// STL includes
#include <vector>

// C/C++ includes
#include <math.h>

// local includes
#include "Det.h"

/** \file Detector.h
 *
 * \brief Configuration of detector layout.
 *
 */

class Detector;

// defined in Detector.cc
extern Detector * gDetector;

/** \class Detector
 *
 * \brief Configuration of detector layout and its misalignment.
 *
 * All values are given in mm and rad, respectively. The detector consists of
 * layers, rods and modules. The default configuration is resembling the TOB
 * Cosmic Rack.
 *
 */
class Detector
{
 public:
  /// default constructor reads detector layout from file
  /// and creates and misaligns all detectors.
  Detector(const char * fname);
  virtual ~Detector();

  /// return number of layers of detector (along z)
  const int    GetNLayer() const { return fNLayer; };
  /// return number of rods in each layer (along x)
  const int    GetNRod() const { return fNRod; };
  /// return number of modules in each rod (along y)
  const int    GetNModule() const { return fNModule; };
  /// distance between layers in z
  const double GetLayerDistance() const { return fLayerDistance; };

  // detector size
  const double GetWidth() const { return fNRod*fModuleWidth; };
  const double GetLength() const { return fNModule*fModuleLength; };
  const double GetHeight() const { return fNLayer*fLayerDistance; };

  /// A misalignment (movement) of the detector. It is given in mm and applied
  /// to either direction (i.e. on + and - direction), being randomized
  /// uniformly.
  const double GetMisalignPosition() const { return fMisalignPosition; };
  /// A misalignment (rotation) of the detector. The angle is given in rad and
  /// applied to either direction, being randomized uniformly.
  const double GetMisalignRotation() const { return fMisalignRotation; };
  /// Get true misalignment in one HepVector
  const CLHEP::HepVector & GetTrueMisalign() const { return fTrueMisalignment; }

  /// Make a linear (array) number for the specified module.
  int MakeID(int layer, int rod, int module) const;

  /// Get detector at that place
  const Det & GetDet(int layer, int rod, int module) const;
  /// reset Detector to its nominal position
  void ResetDet(int moduleID);

 protected:
  // Detector layout
  int    fNLayer;
  int    fNRod;
  int    fNModule;
  double fLayerDistance;
  bool   fRandomPlacement;

  // Module layout
  double fModuleWidth;
  double fModuleLength;
  double fModuleThickness;
  double fModulePitchX;
  double fModulePitchY;

  // Misalignment
  double fMisalignPosition;
  double fMisalignRotation;

  /// Detector collection
  std::vector<Det *> DetVec;

  /// true misalignment parameters
  CLHEP::HepVector fTrueMisalignment;

 private:
  void ReadConfiguration(const char * fname);
  /// Create all detector modules
  void CreateModules();
  /// Create one single detector with ID at specified place
  Det * CreateDetector(int ID, int layer, int rod, int module);
  /// Misalign detector
  void MisalignDetector(Det * adet);
};

#endif
