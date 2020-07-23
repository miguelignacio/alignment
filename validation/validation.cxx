// Includes
//==========

// ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TTree.h"
#include "TLatex.h"

#include <iostream>
#include <fstream>

using namespace std;

// Utilities (INFO, ...)
#include "Utilities.h"

// Plot macro
#include "plot.h"

// linear chi2 fit
#include "linear_least_squares.h"

// global variables, ugly solution to avoid too many parameters...
TH1F * hChi2_ORCA;
TH1F * hChi2_before;
TH1F * hChi2_after;
TH1F * hChi2_other;
TH1F * hChi2_ORCA_per_ndof;
TH1F * hChi2_before_per_ndof;
TH1F * hChi2_after_per_ndof;
TH1F * hChi2_other_per_ndof;
TH1F * hNdof;

// compute chi^2 of tracks, re-compute with given alignment parameters
// this way one can compare chi2 distribution before and after alignment
void compute_chi2(const char * alignmentDataFile, 
		  const TVectorD & alignmentParameters, 
		  int nMax)
{
  // Get information what needs to be aligned
  AlignInfo * info = (AlignInfo *) get_object("AlignInfo", alignmentDataFile);
  if (info == 0) {
    return;
  }
  const int nAlignables = info->GetNAlignables();
  const int nParameters = info->GetNParameters();
  cout << "Computing chi2 for " << nAlignables << " alignables" 
       << " with " << nParameters << " parameters." << endl;

  // sanity check whether this corresponds to input
  if (alignmentParameters.GetNoElements() != nAlignables*nParameters) {
    cout << "Wrong alignmentParameters given, expecte size " << nAlignables*nParameters << endl;
    return;
  }

  // open alignment file and get tree
  TTree * t = init_align_tree(alignmentDataFile);
  AlignEvent * alignEvent = new AlignEvent;
  t->SetBranchAddress("AlignEvent", & alignEvent);

  // loop over entries in tree and compute for each track the chi2
  for (int ii = 0; ii < t->GetEntriesFast(); ii++) {
    if ((nMax > 0) && (ii >= nMax))
      break;
    // read data from file
    t->GetEntry(ii);

    //////////////////////////////////////////////////////////////////////
    // get all data from file
    TVectorD & measurement = *alignEvent->GetMeasurements();
    TArrayI & alignables = *alignEvent->GetIndex();
    const int nTheseAlignables = alignables.GetSize();

    TMatrixDSym & measuredCovariance = *alignEvent->GetMeasuredCovariance();
    TMatrixD & trackDerivatives = *alignEvent->GetTrackDerivatives();
    TVectorD & trackPrediction = *alignEvent->GetTrackPrediction();
    TMatrixD & alignmentDerivatives = *alignEvent->GetAlignmentDerivatives();

    //////////////////////////////////////////////////////////////////////
    // debug output
    PMATRIX("Alignment parameters", alignmentParameters);
    PMATRIX("Alignment derivatives", alignmentDerivatives);
    PMATRIX("Measurement", measurement);
    PMATRIX("Measurement Covariance", measuredCovariance);
    PMATRIX("Track prediction", trackPrediction);
    PMATRIX("Track derivatives", trackDerivatives);

    //////////////////////////////////////////////////////////////////////
    // determine number of hit alignables
    int nHitAlignables=0;
    for (int i = 0; i<nTheseAlignables; i++){
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
    // vector of hit alignment parameters
    TVectorD hitAlignmentParameters(nHitAlignables*nParameters);
    for (int i = 0; i < nHitAlignables; i++) {
      for (int j = 0; j < nParameters; j++)
	hitAlignmentParameters[i*nParameters+j] = alignmentParameters[hitAlignables[i]*nParameters+j];
    }

    //////////////////////////////////////////////////////////////////////
    // track refit

    int ndof = measurement.GetNrows()-4; // track has four parameters

    // build residual
    trackPrediction -= measurement;
    PMATRIX("Residual (before first track refit)", trackPrediction);

    TVectorD trackparameter_corrections = linear_least_squares(trackPrediction, 
							       trackDerivatives, 
							       measuredCovariance);
    PMATRIX("Track parameter corrections", trackparameter_corrections);

    // correct residual by fitted track parameters...
    trackPrediction -= trackDerivatives * trackparameter_corrections;

    PMATRIX("Residual (after first track refit)", trackPrediction);

    //////////////////////////////////////////////////////////////////////
    // compute chi2 and fill histograms

    // ORCA chi2
    hChi2_ORCA->Fill(alignEvent->GetChi2());
    hChi2_ORCA_per_ndof->Fill(alignEvent->GetChi2()/ndof);
    
    // re-compute track chi2 myself
    TMatrixD r(trackPrediction.GetNrows(), 1, trackPrediction.GetMatrixArray());
    TMatrixD rt(1, trackPrediction.GetNrows(), trackPrediction.GetMatrixArray());
    TMatrixD mi(measuredCovariance);
    mi.Invert();
    TMatrixD chi2 = rt * mi * r;
    double d_chi2 = chi2[0][0];
    hChi2_before->Fill(d_chi2);
    hChi2_before_per_ndof->Fill(d_chi2/ndof);
    // PMATRIX("Transposed residual", rt);
    // PMATRIX("inverted covariance matrix", mi);
    // PMATRIX("chi2", chi2);

    // compute chi2 after alignment. Beware of sign! alignmentParameters here
    // are corrections and need to add up with a different sign.
    trackPrediction -= alignmentDerivatives * hitAlignmentParameters;

    PMATRIX("Residual (after alignment corrections are applied", trackPrediction);

    //////////////////////////////////////////////////////////////////////
    // refit track again (since other geometry due to given parameters!)

    trackparameter_corrections = linear_least_squares(trackPrediction, 
						      trackDerivatives,
						      measuredCovariance);
    PMATRIX("Track parameter corrections", trackparameter_corrections);

    // correct residual by fitted track parameters...
    trackPrediction -= trackDerivatives * trackparameter_corrections;

    PMATRIX("Residual (after second track refit)", trackPrediction);

    //////////////////////////////////////////////////////////////////////
    // re-compute track chi2 myself, this time with alignment corrections

    TMatrixD r2(trackPrediction.GetNrows(), 1, trackPrediction.GetMatrixArray());
    TMatrixD rt2(1, trackPrediction.GetNrows(), trackPrediction.GetMatrixArray());
    TMatrixD chi22 = rt2 * mi * r2;
    double d_chi22 = chi22[0][0];
    hChi2_after->Fill(d_chi22);
    hChi2_after_per_ndof->Fill(d_chi22/ndof);

    hNdof->Fill(ndof);
  }
}

void get_parameters(const char * alignmentResultFile, 
		    int nAlignables, int nParameters, TVectorD & alignmentParameters)
{
  // read alignment parameters from alignment file
  TFile histofile(alignmentResultFile, "READ");
  if (!histofile.IsOpen()) {
    cout << "ERR: Opening histo file" << endl;
    return;
  }
  TVectorD(nAlignables*nParameters);
  for (int i = 0; i < nParameters; i++) {
    TH1F * hLimits = (TH1F *) histofile.Get(Form("hLimits%d", i));
    if (hLimits == 0) {
      cout << "ERR: could not read histogram " << i << " from file " << endl;
      return;
    }
    if (nAlignables != hLimits->GetNbinsX()) {
      cout << "ERR: histogram size mismatch, wrong number of params?" << endl;
      cout << "nAlignables = " << nAlignables << ", hLimits->GetNbinsX() = " << hLimits->GetNbinsX() << endl;
      return;
    }
    for (int j = 0; j < nAlignables; j++) {
      alignmentParameters[j*nParameters+i] = hLimits->GetBinContent(j+1);
    }
  }
}

void validate_plot(const char * histogramFile)
{
  setopt(gStyle);
  TFile * histofile = TFile::Open(histogramFile);
  if (histofile == 0) {
    cout << "ERR: opening file" << endl;
    return;
  }
  TH1F * hChi2_ORCA   = (TH1F *) histofile->Get("hChi2_ORCA");
  TH1F * hChi2_before = (TH1F *) histofile->Get("hChi2_before");
  TH1F * hChi2_after  = (TH1F *) histofile->Get("hChi2_after");
  TH1F * hChi2_other  = (TH1F *) histofile->Get("hChi2_other");
  TH1F * hChi2_ORCA_per_ndof   = (TH1F *) histofile->Get("hChi2_ORCA_per_ndof");
  TH1F * hChi2_before_per_ndof = (TH1F *) histofile->Get("hChi2_before_per_ndof");
  TH1F * hChi2_after_per_ndof  = (TH1F *) histofile->Get("hChi2_after_per_ndof");
  TH1F * hChi2_other_per_ndof  = (TH1F *) histofile->Get("hChi2_other_per_ndof");
  if (hChi2_other == 0 || hChi2_after == 0 || hChi2_before == 0 || hChi2_ORCA == 0 ||
      hChi2_other_per_ndof == 0 || hChi2_after_per_ndof == 0 || hChi2_before_per_ndof == 0 || hChi2_ORCA_per_ndof == 0
     ) {
    cout << "ERR: getting histo from file" << endl;
    return;
  }

  TCanvas * c1 = new TCanvas("c1", "Kalman Filter Alignment Validation", 0, 0, 800., 800./TMath::Sqrt(2.));
  setopt(c1);
  setopt(hChi2_before);
  setopt(hChi2_after);
  hChi2_before->Rebin(10);
  hChi2_after->Rebin(10);
  hChi2_before->SetLineStyle(kDotted);
  hChi2_before->SetLineColor(kBlack);
  hChi2_before->SetLineWidth(2);  
  hChi2_before->GetXaxis()->SetRangeUser(0, 100);
  hChi2_before->GetXaxis()->SetTitle("#chi^{2}");
  // cout << "hChi2_before->GetMaximum() " << hChi2_before->GetMaximum() << endl;
  // cout << "hChi2_after->GetMaximum() " << hChi2_after->GetMaximum() << endl;
  hChi2_before->GetYaxis()->SetRangeUser(0, 1.1*TMath::Max(hChi2_before->GetMaximum(), hChi2_after->GetMaximum()));
  hChi2_before->GetYaxis()->SetTitle("number of tracks");
  hChi2_before->GetYaxis()->SetTitleOffset(1.2);
  hChi2_before->Draw();
  hChi2_after->SetLineStyle(kSolid);
  hChi2_after->SetLineColor(kBlue);
  hChi2_after->SetLineWidth(2);  
  hChi2_after->GetXaxis()->SetRangeUser(0, 100);
  hChi2_after->Draw("same");
  TLegend * leg = new TLegend(0.4, 0.6, 0.85, 0.9, "Distribution of track #chi^{2}");
  leg->AddEntry(hChi2_before, Form("before alignment (Mean %4.1f)", hChi2_before->GetMean()), "l");
  leg->AddEntry(hChi2_after, Form("after alignment (Mean %4.1f)", hChi2_after->GetMean()), "l");
  leg->Draw();
  c1->Print("validation_2.pdf");

  TCanvas * c2 = new TCanvas("c2", "Kalman Filter Alignment Validation II", 100, 100, 800., 800./TMath::Sqrt(2.));
  hChi2_before_per_ndof->Rebin(1);
  hChi2_after_per_ndof->Rebin(1);
  hChi2_before_per_ndof->SetLineStyle(kDotted);
  hChi2_before_per_ndof->SetLineColor(kBlack);
  hChi2_before_per_ndof->SetLineWidth(2);  
  hChi2_before_per_ndof->GetXaxis()->SetRangeUser(0, 20);
  hChi2_before_per_ndof->GetXaxis()->SetTitle("#chi^{2}/ndof");
  // cout << "hChi2_before_per_ndof->GetMaximum() " << hChi2_before_per_ndof->GetMaximum() << endl;
  // cout << "hChi2_after_per_ndof->GetMaximum() " << hChi2_after_per_ndof->GetMaximum() << endl;
  hChi2_before_per_ndof->GetYaxis()->SetRangeUser(0, 1.1*TMath::Max(hChi2_before_per_ndof->GetMaximum(), hChi2_after_per_ndof->GetMaximum()));
  hChi2_before_per_ndof->GetYaxis()->SetTitle("number of tracks");
  hChi2_before_per_ndof->GetYaxis()->SetTitleOffset(1.2);
  hChi2_before_per_ndof->Draw();
  hChi2_after_per_ndof->SetLineStyle(kSolid);
  hChi2_after_per_ndof->SetLineColor(kBlue);
  hChi2_after_per_ndof->SetLineWidth(2);  
  hChi2_after_per_ndof->GetXaxis()->SetRangeUser(0, 20);
  hChi2_after_per_ndof->Draw("same");
  leg = new TLegend(0.4, 0.6, 0.85, 0.9, "Distribution of track #chi^{2}/ndof");
  leg->AddEntry(hChi2_before, Form("before alignment (Mean %4.1f)", hChi2_before_per_ndof->GetMean()), "l");
  leg->AddEntry(hChi2_after, Form("after alignment (Mean %4.1f)", hChi2_after_per_ndof->GetMean()), "l");
  leg->Draw();
  c2->Print("validation_3.pdf");
}


void validate(const char * configFile)
{
  // read configuration file
  TEnv mEnv(configFile);

  // file names
  const char * alignmentDataFile = mEnv.GetValue("AlignmentDataFile", 
						 "alignment_data.root");
  const char * alignmentResultFile = mEnv.GetValue("AlignmentResultFile", 
						   "alignment_result.root");
  const char * histogramFile = mEnv.GetValue("HistogramFile", 
					     "alignment_validation.root");

  gLogLevel = mEnv.GetValue("LogLevel", 2);

  // number of tracks
  int nMax = mEnv.GetValue("nMaxTracks", 0);

  // get information on alignment parameters and alignables from alignment data
  AlignInfo * info = (AlignInfo *) get_object("AlignInfo", alignmentDataFile);
  if (info == 0) {
    cout << "Error: Could not read AlignInfo object from file " 
	 << alignmentDataFile << endl;
    return;
  }
  const int nAlignables = info->GetNAlignables();
  const int nParameters = info->GetNParameters();
  cout << "Validation request for " << nAlignables << " alignables" 
       << " with " << nParameters << " parameters." << endl;

  TVectorD alignmentParameters(nAlignables*nParameters);
  get_parameters(alignmentResultFile, nAlignables, nParameters, alignmentParameters);

  hChi2_ORCA            = new TH1F("hChi2ORCA", "#chi^{2} from ORCA tracking (with blow-up)",
				   10000, 0, 1000);
  hChi2_before          = new TH1F("hChi2before", "#chi^{2} from tracks, recomputed by KFA",
				   10000, 0, 1000);
  hChi2_after           = new TH1F("hChi2after", "#chi^{2} from alignment with alignment uncertainty",
				   10000, 0, 1000);
  hChi2_other           = new TH1F("hChi2other", "Other #chi^{2} (for debugging)",
				   10000, 0, 1000);
  hChi2_ORCA_per_ndof   = new TH1F("hChi2ORCA_per_ndof", "#chi^{2}/ndof from ORCA tracking (with blow-up)",
				   10000, 0, 1000);
  hChi2_before_per_ndof = new TH1F("hChi2before_per_ndof", "#chi^{2}/ndof from tracks, recomputed by KFA",
				   10000, 0, 1000);
  hChi2_after_per_ndof  = new TH1F("hChi2after_per_ndof", "#chi^{2}/ndof from alignment with alignment uncertainty",
				   10000, 0, 1000);
  hChi2_other_per_ndof  = new TH1F("hChi2other_per_ndof", "Other #chi^{2}/ndof (for debugging)",
				   10000, 0, 1000);
  hNdof                 = new TH1F("hNdof", "Ndof distribution", 20, 0, 20);

  compute_chi2(alignmentDataFile, alignmentParameters, nMax);
  TLatex * l = 0;
  TCanvas * c1 = new TCanvas("c1", "chi^2", 0, 0, 800, 800./TMath::Sqrt(2.));
  c1->Divide(2,2);
  c1->cd(1);
  hChi2_ORCA->Draw();
  l = new TLatex(0.5, 0.85, "ORCA track #chi^{2}");
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.80, Form("Entries: %7.0f", hChi2_ORCA->GetEntries()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.75, Form("Mean: %7.3f", hChi2_ORCA->GetMean()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.70, Form("RMS:  %7.3f", hChi2_ORCA->GetRMS()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.65, Form("Overflow:  %7.3f", hChi2_ORCA->GetBinContent(hChi2_ORCA->GetNbinsX()+1)));
  l->SetNDC();
  l->Draw();
  c1->cd(2);
  hChi2_before->Draw();
  l = new TLatex(0.5, 0.85, "Kalman #chi^{2} before alignment");
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.80, Form("Entries: %7.0f", hChi2_before->GetEntries()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.75, Form("Mean: %7.3f", hChi2_before->GetMean()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.70, Form("RMS:  %7.3f", hChi2_before->GetRMS()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.65, Form("Overflow:  %7.3f", hChi2_before->GetBinContent(hChi2_before->GetNbinsX()+1)));
  l->SetNDC();
  l->Draw();
  c1->cd(3);
  hChi2_after->Draw();
  l = new TLatex(0.5, 0.85, "Kalman #chi^{2} after alignment");
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.80, Form("Entries: %7.0f", hChi2_after->GetEntries()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.75, Form("Mean: %7.3f", hChi2_after->GetMean()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.70, Form("RMS:  %7.3f", hChi2_after->GetRMS()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.65, Form("Overflow:  %7.3f", hChi2_after->GetBinContent(hChi2_after->GetNbinsX()+1)));
  l->SetNDC();
  l->Draw();
  c1->cd(4);
  // hChi2_other->Draw();
  // l = new TLatex(0.5, 0.85, "Othr #chi^{2} (for debugging)");
  // l->SetNDC();
  // l->Draw();
  // l = new TLatex(0.5, 0.80, Form("Entries: %7.0f", hChi2_other->GetEntries()));
  // l->SetNDC();
  // l->Draw();
  // l = new TLatex(0.5, 0.75, Form("Mean: %7.3f", hChi2_other->GetMean()));
  // l->SetNDC();
  // l->Draw();
  // l = new TLatex(0.5, 0.70, Form("RMS:  %7.3f", hChi2_other->GetRMS()));
  // l->SetNDC();
  // l->Draw();
  // l = new TLatex(0.5, 0.65, Form("Overflow:  %7.3f", hChi2_other->GetBinContent(hChi2_other->GetNbinsX()+1)));
  // l->SetNDC();
  // l->Draw();
  hNdof->Draw();
  l = new TLatex(0.5, 0.85, "Ndof");
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.80, Form("Entries: %7.0f", hNdof->GetEntries()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.75, Form("Mean: %7.3f", hNdof->GetMean()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.70, Form("RMS:  %7.3f", hNdof->GetRMS()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.65, Form("Overflow:  %7.3f", hNdof->GetBinContent(hNdof->GetNbinsX()+1)));
  l->SetNDC();
  l->Draw();


  TCanvas * c2 = new TCanvas("c2", "chi^2/ndof", 100, 100, 800, 800./TMath::Sqrt(2.));
  c2->Divide(2,2);
  c2->cd(1);
  hChi2_ORCA_per_ndof->Draw();
  l = new TLatex(0.5, 0.85, "ORCA track #chi^{2}/ndof");
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.80, Form("Entries: %7.0f", hChi2_ORCA_per_ndof->GetEntries()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.75, Form("Mean: %7.3f", hChi2_ORCA_per_ndof->GetMean()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.70, Form("RMS:  %7.3f", hChi2_ORCA_per_ndof->GetRMS()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.65, Form("Overflow:  %7.3f", hChi2_ORCA_per_ndof->GetBinContent(hChi2_ORCA_per_ndof->GetNbinsX()+1)));
  l->SetNDC();
  l->Draw();
  c2->cd(2);
  hChi2_before_per_ndof->Draw();
  l = new TLatex(0.5, 0.85, "Kalman #chi^{2}/ndof before alignment");
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.80, Form("Entries: %7.0f", hChi2_before_per_ndof->GetEntries()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.75, Form("Mean: %7.3f", hChi2_before_per_ndof->GetMean()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.70, Form("RMS:  %7.3f", hChi2_before_per_ndof->GetRMS()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.65, Form("Overflow:  %7.3f", hChi2_before_per_ndof->GetBinContent(hChi2_before_per_ndof->GetNbinsX()+1)));
  l->SetNDC();
  l->Draw();
  c2->cd(3);
  hChi2_after_per_ndof->Draw();
  l = new TLatex(0.5, 0.85, "Kalman #chi^{2}/ndof after alignment");
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.80, Form("Entries: %7.0f", hChi2_after_per_ndof->GetEntries()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.75, Form("Mean: %7.3f", hChi2_after_per_ndof->GetMean()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.70, Form("RMS:  %7.3f", hChi2_after_per_ndof->GetRMS()));
  l->SetNDC();
  l->Draw();
  l = new TLatex(0.5, 0.65, Form("Overflow:  %7.3f", hChi2_after_per_ndof->GetBinContent(hChi2_after_per_ndof->GetNbinsX()+1)));
  l->SetNDC();
  l->Draw();
  c2->cd(4);
  // hChi2_other_per_ndof->Draw();
  // l = new TLatex(0.5, 0.85, "Othr #chi^{2}/ndof (for debugging)");
  // l->SetNDC();
  // l->Draw();
  // l = new TLatex(0.5, 0.80, Form("Entries: %7.0f", hChi2_other_per_ndof->GetEntries()));
  // l->SetNDC();
  // l->Draw();
  // l = new TLatex(0.5, 0.75, Form("Mean: %7.3f", hChi2_other_per_ndof->GetMean()));
  // l->SetNDC();
  // l->Draw();
  // l = new TLatex(0.5, 0.70, Form("RMS:  %7.3f", hChi2_other_per_ndof->GetRMS()));
  // l->SetNDC();
  // l->Draw();
  // l = new TLatex(0.5, 0.65, Form("Overflow:  %7.3f", hChi2_other_per_ndof->GetBinContent(hChi2_other_per_ndof->GetNbinsX()+1)));
  // l->SetNDC();
  // l->Draw();
  // create output histogram file and histograms
  // also store alignable and parameter information in the 
  c1->Print("validation_1.pdf");
  TFile histofile(histogramFile, "RECREATE");
  hChi2_ORCA->Write("hChi2_ORCA");
  hChi2_before->Write("hChi2_before");
  hChi2_after->Write("hChi2_after");
  hChi2_other->Write("hChi2_other");
  hChi2_ORCA_per_ndof->Write("hChi2_ORCA_per_ndof");
  hChi2_before_per_ndof->Write("hChi2_before_per_ndof");
  hChi2_after_per_ndof->Write("hChi2_after_per_ndof");
  hChi2_other_per_ndof->Write("hChi2_other_per_ndof");
  hNdof->Write("hNdof");
  histofile.Close();
  validate_plot(histogramFile);
}

// Entry point
//=============
int main(int argc, char * argv[])
{
  cout << 
    "Kalman Alignment - track based alignment of high-energy particle physics detectors.\n"
    "(validation package)\n"
    "\n"
    "Copyright (C) 2005-2011 Martin Weber\n"
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
	 << "validate ConfigFile" << endl
	 << endl
	 << "ConfigFile: Contains the configuration for validation" << endl;
    exit(1);
  }

  validate(argv[1]);
}
