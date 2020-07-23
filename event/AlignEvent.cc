#include <iostream>

using namespace std;

#include "AlignEvent.h"
#include "TMath.h"

class RecHit;
class DetUnit;

ClassImp(AlignInfo)
ClassImp(AlignEvent)

//////////////////////////////////////////////////////////////////////
// AlignInfo

AlignInfo::AlignInfo()
  : fNAlignables(0), fNParameters(0)
{
}

AlignInfo::AlignInfo(int nAlignables, int nParameters)
  : fNAlignables(nAlignables),
    fNParameters(nParameters)
{ 
}

AlignInfo::~AlignInfo()
{
}

//////////////////////////////////////////////////////////////////////
// AlignEvent

AlignEvent::AlignEvent()
{
  fIndex = new TArrayI;
  fMeasurements = new TVectorD;
  fMeasuredCovariance = new TMatrixDSym;
  fTrackDerivatives = new TMatrixD;
  fAlignmentDerivatives = new TMatrixD;
  fTrackPrediction = new TVectorD;
}

AlignEvent::~AlignEvent()
{
  delete fIndex;
  delete fMeasurements;
  delete fMeasuredCovariance;
  delete fTrackDerivatives;
  delete fAlignmentDerivatives;
  delete fTrackPrediction;
}
