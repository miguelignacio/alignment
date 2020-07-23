// Includes
//==========

// ROOT
#include "TArrayD.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TEnv.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRint.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVectorD.h"

#include <iostream>
#include <fstream>

// Utilities (INFO, ...)
#include "Utilities.h"
#include "AlignEvent.h"

using namespace std;

// TODO list
//===========
// - multi-track events?
// - symmetry bug?
// - optimise loops


// New Kalman Filter Alignment Algorithm
//=======================================
void kfa(const char * configFile = "align.cfg", int gLogLevel = 2)
{
  int fileNr = 0;

  // read config file
  //------------------
  TEnv mEnv(configFile);

  // number of tracks
  int nMax = mEnv.GetValue("nMaxTracks", 0);

  // file names
  const char * inputFileName = mEnv.GetValue("inputFileName", "alignment_data.root");
  const char * outputFileName = mEnv.GetValue("outputFileName", "alignment_result.root");
  const char * outlierRejectionFileName = mEnv.GetValue("outlierRejectionFileName", "outliers.log");

  // deterministic annealing
  double AnnealingFactor = mEnv.GetValue("annealingFactor", 1000.);
  int    AnnealingEvents = mEnv.GetValue("annealingEvents", 0);  

  // alignable fixation
  // number of fixed alignables
  int    nFixedAlignables = 0;
  while (mEnv.Defined(Form("fixedAlignable%d", nFixedAlignables + 1))){
    nFixedAlignables++;
  }
  int    fixedAlignable[nFixedAlignables];
  for (int i=0; i < nFixedAlignables; i++) {
    fixedAlignable[i] = mEnv.GetValue(Form("fixedAlignable%d", i+1), 1);
  }

  // parameter fixation
  // number of fixed parameters
  int    nFixedParameters = 0;
  while (mEnv.Defined(Form("fixedParameter%d", nFixedParameters + 1))){
    nFixedParameters++;
  }
  int    fixedParameter[nFixedParameters];
  for (int i=0; i < nFixedParameters; i++) {
    fixedParameter[i] = mEnv.GetValue(Form("fixedParameter%d", i+1), 1);
  }

  // fixation error
  double fixationError = mEnv.GetValue("fixationError", 1E-10);

  // start value
  double startValue = mEnv.GetValue("startValue", 0.);
  double startError = mEnv.GetValue("startError", 1E-2);
  double startFixedValue = mEnv.GetValue("startFixedValue", 0.);

  // outlier rejection
  double probCut = mEnv.GetValue("probabilityCut", 0.5);
  double deviationCut = mEnv.GetValue("deviationCut", 1.0);

  // random order
  bool   randomOrder = mEnv.GetValue("randomEventOrder", true);
  int    randomOrderSeed = mEnv.GetValue("randomEventOrderSeed", 8071984);


  // Get information what needs to be aligned
  //------------------------------------------
  AlignInfo * info = (AlignInfo *) get_object("AlignInfo", inputFileName);

  if (info == 0) {
    cout << "Error: Could not read AlignInfo object." << endl;
    return;
  }
  const int nAlignables = info->GetNAlignables();
  const int nParameters = info->GetNParameters();
  cout << "Alignment request for " << nAlignables << " alignables" 
       << " with " << nParameters << " parameters." << endl;

  // init number of updates counters for each alignable
  int nUpdates[nAlignables];
  for (int i=0; i < nAlignables; i++){
    nUpdates[i]=0;
  }

  // open alignment file and get tree
  TTree * t = init_align_tree(inputFileName);
  AlignEvent * alignEvent = new AlignEvent;
  t->SetBranchAddress("AlignEvent", & alignEvent);

  // number of events to be processed (-1 or 0 means all)
  if (nMax <= 0)
    nMax=t->GetEntriesFast()-1;

  cout << "Processing " << (nMax + 1) << " events" << endl;

  // create output histogram file and histograms
  // also store alignable and parameter information in the file
  // TODO: configure file path via cfg file
  TFile histofile(outputFileName, "RECREATE");
  TH1F * hAliPar[nAlignables*nParameters];
  for (int i = 0; i < nAlignables*nParameters; i++) {
    hAliPar[i] = new TH1F(Form("hAliPar%d", i), 
			  Form("Evolution of alignment parameter %d", i),
			  nMax, 0, nMax);
  }
  TH1F * hChi2_Track = new TH1F("hChi2Track", "#chi^{2} from tracking",
				1000, 0, 1000);
  TH1F * hChi2_Ali   = new TH1F("hChi2Ali", "#chi^{2} from alignment",
				10000, 0, 10000);
  TH1F * hChi2_Ali2  = new TH1F("hChi2Ali2", "#chi^{2} from alignment with alignment uncertainty",
				1000, 0, 1000);

  /*
  // store nAlignables, nParameters
  // cannot store just an integer, so use TNamed instead
  TNamed * nameNAlignables = new TNamed("nAlignables", Form("%d", nAlignables));
  histofile.WriteObject(nameNAlignables, "nAlignables", "WriteDelete");
  TNamed * nameNParameters = new TNamed("nParameters", Form("%d", nParameters));
  histofile.WriteObject(nameNParameters, "nParameters", "WriteDelete");
  */

  // Create alignment store
  TVectorD    alignmentParameters(nAlignables*nParameters);
  TVectorD    copyAlignmentParameters(nAlignables*nParameters);


  // Print info about alignment options
  //------------------------------------
  cout << "Track probability cut: " << probCut << endl;
  if (deviationCut > 0.0) { cout << "Deviation cut: " << deviationCut << endl; }
  if (AnnealingEvents != 0) { cout << "Deterministic annealing until event " << AnnealingEvents << endl; }
  if (mEnv.GetValue("chi2Ali2Cut", -1) != -1) { cout << "Chi2Ali2Cut: " << mEnv.GetValue("chi2Ali2Cut", -1) << endl; }

  int nChi2Ali2CutEvents = 0; // counter for cut events


  // Initialize parameters and covariance
  //--------------------------------------
  cout << "Start value: " << startValue << endl;
  for (int i = 0; i < nAlignables*nParameters; i++) {
    alignmentParameters[i] = startValue;
    // could also use random starting values
    //alignmentParameters[i] = gRandom->Uniform(-startval, startval);
  }
  cout << "Start error: " << startError << endl;
  TMatrixDSym alignmentCovariance(nAlignables*nParameters);
  TMatrixDSym copyAlignmentCovariance(nAlignables*nParameters);
  for (int i = 0; i < nAlignables*nParameters; i++) {
    alignmentCovariance[i][i] = startError;
  }

  // Fix alignables
  //----------------
  for (int i=0; i<nFixedAlignables; i++){
    if (fixedAlignable[i] <= nAlignables){
      if (gLogLevel>=2)
	cout << "Fixing alignable " << fixedAlignable[i] << " to " << startFixedValue
	     << " +- " << fixationError << endl;

      for (int j=0; j<nParameters; j++){
	alignmentParameters[(fixedAlignable[i]-1)*nParameters + j] = startFixedValue;
	alignmentCovariance[(fixedAlignable[i]-1)*nParameters + j][(fixedAlignable[i]-1)*nParameters + j] = fixationError;
      }
    }else{
      if (gLogLevel>=2)
	cout << "Warning: Could not fix alignable " << fixedAlignable[i] << endl;
    }
  }

  // Fix parameters
  //----------------
  for (int i=0; i<nFixedParameters; i++){
    if (fixedParameter[i] < nParameters){
      if (gLogLevel>=2)
	cout << "Fixing parameter " << fixedParameter[i] << " to 0"
	     << " +- " << fixationError << endl;

      for (int j=0; j<nAlignables; j++){
	alignmentParameters[j*nParameters + fixedParameter[i]] = 0;
	alignmentCovariance[j*nParameters + fixedParameter[i]][j*nParameters + fixedParameter[i]] = fixationError;
      }
    }else{
      if (gLogLevel>=2)
	cout << "Warning: Could not fix parameter " << fixedParameter[i] << endl;
    }
  }

  // Fix specific alignable in a specific parameter
  //------------------------------------------------
  for (int iAlignable=0; iAlignable < nAlignables; iAlignable++){
    for (int iParameter=0; iParameter < nParameters; iParameter++){
      if (mEnv.GetValue(Form("fixParameter%dofAlignable%d", iParameter, iAlignable + 1), false)){
	if (gLogLevel>=2)
	  cout << "Fixing parameter " << iParameter << " of alignable " << iAlignable + 1 << " to " << startFixedValue
	       << " +- " << fixationError << endl;

	  alignmentParameters[iAlignable*nParameters + iParameter] = startFixedValue;
	  alignmentCovariance[iAlignable*nParameters + iParameter][iAlignable*nParameters + iParameter] = fixationError;
      }
    }
  }

  // Create random map for track loop
  //----------------------------------
  int iMap[t->GetEntriesFast()];
  for (int ii = 0; ii < t->GetEntriesFast(); ii++){
    iMap[ii]=ii;
  }
  if (randomOrder){
    cout << "Using random event order with seed " << randomOrderSeed << endl;
    TRandom rnd;
    rnd.SetSeed(randomOrderSeed);
    for (int ii=0; ii<t->GetEntriesFast(); ii++){
      int temp=iMap[ii];
      int iRnd=(int)rnd.Uniform(0,t->GetEntriesFast());
      iMap[ii]=iMap[iRnd];
      iMap[iRnd]=temp;
    }
  }
  
  // outlier rejection log output
  ofstream logfileOutliers(outlierRejectionFileName);

  // counters
  int nStep = 0;
  int nOutliers1 = 0;
  int nOutliers2 = 0;


  // Loop
  //------
  for (int ii = 0; ii < t->GetEntriesFast(); ii++) {
    // check max events
    if (nStep >= nMax)
      break;

    t->GetEntry(iMap[ii]);

    // Outlier rejection II
    //----------------------
    // save alignment parameters and alignment covariance
    if (deviationCut > 0.0){
      copyAlignmentParameters = alignmentParameters;
      copyAlignmentCovariance = alignmentCovariance;
    }

    // get all data from file
    //------------------------
    TVectorD & measurement = *alignEvent->GetMeasurements();
    const int nMeasurements = measurement.GetNrows();
    TArrayI & alignables = *alignEvent->GetIndex();
    const int nTheseAlignables = alignables.GetSize();

    TMatrixDSym & measuredCovariance = *alignEvent->GetMeasuredCovariance();
    TMatrixD & trackDerivatives = *alignEvent->GetTrackDerivatives();

    TMatrixD & alignmentDerivatives = *alignEvent->GetAlignmentDerivatives();

    // Skip track if there are more than one hits in the same alignable
    // TODO: prepare algorithm to work on these tracks!
//     if (nParameters*nMeasurements/2!=alignmentDerivatives.GetNcols())
//       continue;

    TVectorD & trackPrediction = *alignEvent->GetTrackPrediction();


    // Some generic information about event before quality cuts are applied
    //----------------------------------------------------------------------
    INFO("Run " << alignEvent->GetRun() << " event " << alignEvent->GetEvent());
    INFO("Step number " << nStep);
    INFO("Tree entry " << ii);
    INFO("Chi2  " << alignEvent->GetChi2());
    INFO("Ndof  " << alignEvent->GetNdof());
    double prob = TMath::Prob(alignEvent->GetChi2(), alignEvent->GetNdof());
    INFO("Prob  " << prob);
    hChi2_Track->Fill(alignEvent->GetChi2());

    // Data quality checks
    //---------------------
    // count stereo measurements (by covariance diagonal entries)
    int nStereo = 0;
    for (int i = 0; i < nMeasurements/2; i++) {
      if (measuredCovariance[i*2+1][i*2+1] < 5)
	nStereo++;
    }
    INFO("Found " << nStereo << " stereo measurements");
    // if (skipStereoModules && nStereo > 0) {
    //   INFO("Skipping event because of stereo modules");
    //   continue;
    // }
    if (prob < probCut) { 
      INFO("Skipping event because of track probability");
      nOutliers1++;
      continue;
    }

    // TODO: optimize the following loops
    // determine number of hit alignables
    int nHitAlignables=0;
    for (int i = 0; i<nTheseAlignables; i++){
      // update update counter
      nUpdates[alignables[i]]++;

      bool haveAlready=false;
      for (int j = 0; j<i; j++){
	if (alignables[j]==alignables[i])
	  haveAlready=true;	
      }

      if (!haveAlready)
	nHitAlignables++;
    }
    
    // create vector of hit alignables
    int hitAlignables[nHitAlignables];    
    int counter=0;
    for (int i = 0; i<nTheseAlignables; i++){
      bool haveAlready=false;
      for (int j = 0; j<i; j++){
	if (alignables[j]==alignables[i])
	  haveAlready=true;	
      }

      if (!haveAlready){
	hitAlignables[counter]=alignables[i];
	counter++;	
      }
    }    

    // Debug output
    //--------------
    if (gLogLevel > 2) { 
      cout << "Alignable IDs" << endl;
      for (int i = 0; i < nTheseAlignables; i++) {
	cout << alignables[i] << endl;
      }
      cout << nHitAlignables << " alignables hit" << endl;
    }
    PMATRIX("Alignment parameters", alignmentParameters);
    PMATRIX("Alignment covariance", alignmentCovariance);
    PMATRIX("Alignment derivatives", alignmentDerivatives);
    PMATRIX("Measurement", measurement);
    PMATRIX("Measurement Covariance", measuredCovariance);
    PMATRIX("Track prediction", trackPrediction);
    PMATRIX("Track derivatives", trackDerivatives);

    // Deterministic Annealing (geometric cooling scheme)
    //----------------------------------------------------
    if (nStep < AnnealingEvents) {
      INFO("Annealing factor: " << TMath::Power(AnnealingFactor, (AnnealingEvents-nStep)/(static_cast<float>(AnnealingEvents))));
      for (int i=0; i < nMeasurements/2; i++){
	bool isFixed = false;

	// check if measurement is on a fixed alignable
	for (int iFixed=0; iFixed < nFixedAlignables; iFixed++){
	  if (alignables[i] == fixedAlignable[iFixed])
	    isFixed=true;
	}

	// do not apply annealing if measurement is on a fixed alignable
	if (!isFixed){
	  measuredCovariance[2*i][2*i] *= TMath::Power(AnnealingFactor, (AnnealingEvents-nStep)/(static_cast<float>(AnnealingEvents)));
	  measuredCovariance[2*i + 1][2*i + 1] *= TMath::Power(AnnealingFactor, (AnnealingEvents-nStep)/(static_cast<float>(AnnealingEvents)));
	}
      }
      PMATRIX("Annealed covariance ", measuredCovariance);
    }
    else {
      static bool said = false;
      if (!said) {
 	INFO("Annealing stopped - max number of events reached");
	said = true;
      }
    }


    // vector of hit alignables
    TVectorD hitAlignmentParameters(nHitAlignables*nParameters);
    for (int i = 0; i < nHitAlignables; i++) {
      for (int j = 0; j < nParameters; j++)
	hitAlignmentParameters[i*nParameters+j] = alignmentParameters[hitAlignables[i]*nParameters+j];
    }

    // compute residual (m - c - Ad) in Fruehwirth's notation
    trackPrediction -= measurement + alignmentDerivatives * hitAlignmentParameters;

    // Compute chi2 myself... track only
    //-----------------------------------
    PMATRIX("Residual", trackPrediction);
    TMatrixD r(trackPrediction.GetNrows(), 1, trackPrediction.GetMatrixArray());
    PMATRIX("r ", r);
    TMatrixD rt(1, trackPrediction.GetNrows(), trackPrediction.GetMatrixArray());
    PMATRIX("rt ", rt);
    TMatrixD mi(measuredCovariance);
    mi.Invert();
    TMatrixD chi2 = rt * mi * r;
    double d_chi2 = chi2[0][0];
    INFO("TTrack chi2 = " << alignEvent->GetChi2() << " my chi2 (track only) = " << d_chi2);
    hChi2_Ali->Fill(d_chi2);

    // Kalman filter alignment
    //-------------------------
    // we need the matrix D*A^T in the Kalman Filter formalism as a temporary
    // matrix storage. This matrix is used in both the update of the
    // parameters and of the covariance matrix. It is the largest matrix we
    // use in KFA, so it might be worth to improve computation speed here.
    TMatrixD DAT(nAlignables*nParameters, nMeasurements);
    
    // Matrix multiplication
    for(int i=0; i<nAlignables*nParameters; i++){
      for (int j=0; j<nMeasurements; j++){
	// multiplication for matrix element [i;j] (line i, column j)
	double val=0;

	// skip parameters not affected by current track.
	// alignmentDerivatives only contains affected parameters. Their ID can be
	// looked up in alignables. Determine their column number in alignmentCovariance
	// that way.
	for (int k=0; k< nHitAlignables; k++){
	  for (int l=0; l<nParameters; l++){
	    val+=alignmentCovariance[i][hitAlignables[k]*nParameters+l]*alignmentDerivatives[j][k*nParameters+l];
	  }
	}
	DAT[i][j]=val;	
      }
    }
    PMATRIX("DAT = ", DAT);


    // now we want to compute V_D and (V_D)^-1 (which I call "VDI"), in order
    // to compute the Sherman-Morrison inverse and the matrix G (which I call
    // "Weight") of the filter.

    // initialize matrix with measurement covariance
    TMatrixDSym VDI(measuredCovariance);
    // add by hand the multiplication with alignmentDerivatives since
    // ROOT does not allow me to do the multiplication (end product is symmetric)
    // Improved multiplication
    for (int line = 0; line < nMeasurements; line++) {
      for (int column = line; column < nMeasurements; column++) {
	// build matrix product ADAT
	double val = 0;
	for (int kAlignable = 0; kAlignable < nHitAlignables; kAlignable++) {
	  for (int kParameter = 0; kParameter < nParameters; kParameter++){
	    int k1=kAlignable*nParameters + kParameter; // unsorted index
	    int k2=hitAlignables[kAlignable]*nParameters + kParameter; // sorted index
	    val += alignmentDerivatives[line][k1]*DAT[k2][column];
	  }
	}
	VDI[line][column] += val;
	/// \bug : Why can't I use just the symmetry? Seems to be a ROOT
	/// problem...
	VDI[column][line] = VDI[line][column];
      }
    }

    PMATRIX("VD", VDI);
    
    // calculate eigenvalues
    if (gLogLevel > 2){
      TVectorD vecEValues;
      PMATRIX("VD Eigenvectors: ", VDI.EigenVectors(vecEValues));
      PMATRIX("VD Eigenvalues: ", vecEValues);
    }

    VDI.Invert();
    PMATRIX("VDI", VDI);

    TMatrixDSym Weight(VDI);

    // compute chi^2 including the uncertainty from alignment and track
    chi2 = rt * VDI * r;
    d_chi2 = chi2[0][0];
    hChi2_Ali2->Fill(d_chi2);
    INFO("TTrack chi2 = " << alignEvent->GetChi2() << " my chi2 with ali = " << d_chi2);

    Weight.SimilarityT(trackDerivatives);
    PMATRIX("Weight", Weight);

    // calculate eigenvalues
    if (gLogLevel > 2){
      TVectorD vecEValues;
      PMATRIX("Weight Eigenvectors: ", Weight.EigenVectors(vecEValues));
      PMATRIX("Weight Eigenvalues: ", vecEValues);
    }

    Weight.Invert();
    Weight.Similarity(trackDerivatives);
    Weight.Similarity(VDI);
    Weight *= -1;
    Weight += VDI;
    PMATRIX("Weight matrix after Sherman Morrison inverse", Weight);

    // full update of alignment parameters
    // the next line is neither CPU nor memory effective but works.
    alignmentParameters += DAT * Weight * trackPrediction;
    PMATRIX("Updated alignment parameters:", alignmentParameters);

    // full update of covariance matrix
    // the next line is neither CPU nor memory effective but works.
    Weight.Similarity(DAT);
    alignmentCovariance -= Weight;
    PMATRIX("Updated alignment covariance", alignmentCovariance);

    // Outlier rejection II
    //----------------------
    // check for exceeding deviation and eventually restore alignment parameters and covariance
    if (deviationCut > 0.0){
      bool filterEvent = false;
      for (int i=0; i < nAlignables*nParameters; i++){
	if (TMath::Abs(alignmentParameters[i] - copyAlignmentParameters[i]) > deviationCut * TMath::Sqrt(copyAlignmentCovariance[i][i])){
	  //cout << i << ": " << TMath::Abs(alignmentParameters[i] - copyAlignmentParameters[i]) << ", " << TMath::Sqrt(copyAlignmentCovariance[i][i]) << endl;
	  filterEvent = true;
	}
      }
      if (filterEvent){
	// restore
	INFO("Outlier rejection stage 2: Skipping event, restoring parameters and covariance");
	alignmentParameters = copyAlignmentParameters;
	alignmentCovariance = copyAlignmentCovariance;
	nOutliers2++;

	logfileOutliers << "Filtering: run " << alignEvent->GetRun()
			<< ", event " << alignEvent->GetEvent()
			<< ", step " << nStep
			<< ", chi2 " << alignEvent->GetChi2()
			<< ", Ndof " << alignEvent->GetNdof()
			<< "\n";
      }
    }

    // fill histograms
    for (int i = 0; i < nAlignables*nParameters; i++) {
      hAliPar[i]->SetBinContent(nStep+1, alignmentParameters[i]);
      hAliPar[i]->SetBinError(nStep+1, TMath::Sqrt(alignmentCovariance[i][i]));
    }

    // increase alignment step counter
    //---------------------------------
    nStep++;
    // occasionally give status report
    if (nStep % 10000 == 0){
      TDatime theTime = TDatime();
      cout << theTime.GetHour() << ":" << theTime.GetMinute() << "'" << theTime.GetSecond()
	   << ": " << nStep << " tracks finished ..." << endl;
    }
  }

  // Store limits histograms
  //-------------------------
  for (int i=0; i<nParameters; i++){
    TH1F * hLimits = new TH1F(Form("hLimits%d", i), Form("Evolution limits of parameter %d", i),
			      nAlignables, 0, nAlignables);

    for (int j=0; j<nAlignables; j++){
      hLimits->SetBinContent(j+1, alignmentParameters[j*nParameters + i]);
      hLimits->SetBinError(j+1, TMath::Sqrt(alignmentCovariance[j*nParameters + i][j*nParameters + i]));
      if (gLogLevel >= 3){
	// output of parameter limits
	cout << "align. " << j << ", par. " << i << ": " <<  alignmentParameters[j*nParameters + i]
	     << " +- " << TMath::Sqrt(alignmentCovariance[j*nParameters + i][j*nParameters + i]) << endl;
      }
    }
  }

  alignmentCovariance.Write("CovarianceMatrix");


  // Store number of updates histograms
  //------------------------------------
  TH1F * hNOfUpdates = new TH1F("hNOfUpdates", "Hits per Alignable", nAlignables, 0.5, 0.5 + nAlignables);
  for (int j=0; j<nAlignables; j++){
    hNOfUpdates->SetBinContent(j+1, nUpdates[j]);
  }

  histofile.Write();
  histofile.Close();

  cout << nOutliers1 << " events rejected by outlier rejection stage 1" << endl;
  cout << nStep << " tracks processed" << endl;
  if (deviationCut > 0.0){
    cout << nOutliers2 << " events ignored by outlier rejection stage 2" << endl;
  }

  logfileOutliers.close();
}


// Entry point
//=============
int main(int argc, char * argv[])
{
  cout << 
    "Kalman Alignment - track based alignment of high-energy particle physics detectors.\n"
    "(alignment package)\n"
    "\n"
    "Copyright (C) 2005-2011 Martin Weber, Daniel Sprenger\n"
    "\n"
    "This program is free software: you can redistribute it and/or modify\n"
    "it under the terms of the GNU General Public License as published by\n"
    "the Free Software Foundation, either version 3 of the License, or\n"
    "(at your option) any later version.\n"
    "\n"
    "This program is distributed in the hope that it will be useful,\n"
    "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
    "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
    "GNU General Public License for more details.\n"
    "\n"
    "You should have received a copy of the GNU General Public License\n"
    "along with this program.  If not, see <http://www.gnu.org/licenses/>.\n";

  // check command line arguments
  if (argc != 2) {
    cout << "Usage:" << endl
	 << endl
	 << "align ConfigFile" << endl
	 << endl
	 << "ConfigFile: Contains the configuration for alignment parameters," << endl
	 << "            verbosity and more." << endl
	 << endl;
    exit(1);
  }

  kfa(argv[1]);
}

