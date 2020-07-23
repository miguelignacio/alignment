#ifndef _AlignEvent_h
#define _AlignEvent_h

#include "TObject.h"
#include "TClonesArray.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TArrayI.h"

class AlignInfo : public TObject {
  int fNAlignables; // number of alignables
  int fNParameters; // number of alignment parameters
  
 public:
  AlignInfo();
  AlignInfo(int nAlignables, int nParameters);
  int GetNAlignables() const { return fNAlignables; };
  int GetNParameters() const { return fNParameters; };
  virtual ~AlignInfo();  
  
  ClassDef(AlignInfo, 1); // Alignment Information
};

class AlignEvent {
protected:
  int           fRun;          // run number
  int           fEvent;        // event number
  double        fChi2;         // chi^2 of track fit
  int           fNdof;         // number of degrees of freedom
  int           fTrackNumber;  // number of corresponding track
  TArrayI     * fIndex;                //->index of Alignables
  TVectorD    * fMeasurements;         //->measurements
  TMatrixDSym * fMeasuredCovariance;   //->covariance of measurements
  TMatrixD    * fTrackDerivatives;     //->track derivatives
  TMatrixD    * fAlignmentDerivatives; //->alignment derivatives
  TVectorD    * fTrackPrediction;      //->track prediction

 public:
  AlignEvent();
  virtual ~AlignEvent();
  void SetRun(int run) { fRun = run; };
  void SetEvent(int event) { fEvent = event; };
  void SetChi2(double chi2) { fChi2 = chi2; };
  void SetNdof(int ndof) { fNdof = ndof; };
  void SetTrackNumber(int number) { fTrackNumber = number; };
  int  GetRun() const { return fRun; };
  int  GetEvent() const { return fEvent; };
  double GetChi2() const { return fChi2; };
  int  GetNdof() const { return fNdof; };
  int  GetTrackNumber() const { return fTrackNumber; };
  // const getters for index, vectors, matrices
  const TArrayI     * GetIndex() const { return fIndex; };
  const TVectorD    * GetMeasurements() const { return fMeasurements; };
  const TMatrixDSym * GetMeasuredCovariance() const { return fMeasuredCovariance; };
  const TMatrixD    * GetTrackDerivatives() const { return fTrackDerivatives; };
  const TMatrixD    * GetAlignmentDerivatives() const { return fAlignmentDerivatives; };
  const TVectorD    * GetTrackPrediction() const { return fTrackPrediction; };

  // non-const getters for index, vectors, matrices to allow setting
  // without constructing
  TArrayI     * GetIndex() { return fIndex; };
  TVectorD    * GetMeasurements() { return fMeasurements; };
  TMatrixDSym * GetMeasuredCovariance() { return fMeasuredCovariance; };
  TMatrixD    * GetTrackDerivatives() { return fTrackDerivatives; };
  TMatrixD    * GetAlignmentDerivatives() { return fAlignmentDerivatives; };
  TVectorD    * GetTrackPrediction() { return fTrackPrediction; };

  ClassDef(AlignEvent, 1); // Alignment event
};

#endif
