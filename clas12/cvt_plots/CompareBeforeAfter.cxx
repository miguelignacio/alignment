#include <cstdlib>
#include <list>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TStyle.h>
#include <TLine.h>
#include <TText.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <TLegend.h>
#include <TPaveText.h>
#include "../../event/AlignEvent.h"


#include <iostream>



#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"


using namespace std;

double dmin=-30,dmax = 30;
double phimin=-TMath::Pi(), phimax = TMath::Pi();
double tandipmin= -1,tandipmax =1.5;
double zmin=-70, zmax =10;
TString plotsDir ="";

TF1* g = new TF1("gaus","gaus", -.1,.1);
double getSigma(TH1* h){
  double mu = h->GetMean();
  double sigma = h->GetStdDev();
  
  double minfit = mu-2*sigma;
  double maxfit = mu+2*sigma;
  g->SetParameter("Mean", mu);
  g->SetParameter("Sigma", sigma);
  h->Fit(g, "N", "", minfit, maxfit);
  mu = g->GetParameter("Mean");
  sigma = g->GetParameter("Sigma");
  return sigma;
}
TGraphErrors* createProfile(TH2D * h, int color, int markerstyle, double shift,double ywindow=1.5){
  
  int n =h->GetXaxis()->GetNbins();
  double x[n];
  double y[n];
  double ex[n];
  double ey[n];

  double largestMu = 0;
  double largestMuError = 0;
  double chi2 = 0;
  for(int i = 0; i<n;i++){
    TH1* proj = h->ProjectionY("temp", i,i+1);
    double mu = proj->GetMean();
    //double dmu = proj->GetMeanError();
    double sigma = proj->GetStdDev();
    
    double minfit = mu-2*sigma;
    double maxfit = mu+2*sigma;
    g->SetParameter("Mean", mu);
    g->SetParameter("Sigma", sigma);
    proj->Fit(g, "NQ", "", minfit, maxfit);
    mu = g->GetParameter("Mean");
    sigma = g->GetParameter("Sigma");
    double dmu = g->GetParError(1);
    if(sigma>5){
      mu = proj->GetMean();
      sigma = proj->GetStdDev();
      //proj->Draw();
      //gPad->SaveAs(plotsDir+ "/error.pdf");
      
    }
    x[i] = h->GetXaxis()->GetBinCenter(i+1)+shift;
    ex[i] = 0;
    y[i] = mu;
    ey[i] = sigma;
    if(abs(mu)>largestMu){
      largestMu = abs(mu);
      largestMuError = dmu;
      
    }
    //for empty or nearly empty bins, hide the marker
    if(proj->GetEntries()<10){
      y[i]=100;
    }
    chi2+=mu*mu/(dmu*dmu);
  }
  
  TGraphErrors * profile = new TGraphErrors(n, x, y, ex, ey);
  profile->SetMarkerStyle(markerstyle);
//profile->SetMarkerSize(1);
  profile->SetMarkerColor(color);
  profile->SetLineColor(color);
  profile->SetLineStyle(0);
  profile->GetXaxis()->SetRangeUser(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  profile->SetTitle(h->GetTitle());
  profile->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  profile->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
  profile->GetYaxis()->SetTitleOffset(h->GetYaxis()->GetTitleOffset());
  
  profile->GetHistogram()->SetMaximum(ywindow);   // along
  profile->GetHistogram()->SetMinimum(-ywindow);  //   Y
  
  cout << h->GetName() << " worst mu: " << largestMu << "+-"<< largestMuError << endl;
  cout << h->GetName() << " chi2: " << chi2 << endl;
  return profile;
//return h->ProfileX(h->GetName()+(TString&)"_prof");
  
}

int main(int argc, char * argv[]) {
  // Record start time
  //auto start = std::chrono::high_resolution_clock::now();
  
  TString inputFileNameBefore;
  TString inputFileNameAfter;
  //TString plotsDir;

  double maxResid=5;
  vector<TString> inputFiles;
  TString label ="";
  bool isMC = 0;
  for(Int_t i=1;i<argc;i++){
    TString opt=argv[i];
    if((opt.Contains("--before="))){
      inputFileNameBefore=opt(9,opt.Sizeof());
    }
    if((opt.Contains("--after="))){
      inputFileNameAfter=opt(8,opt.Sizeof());
    }
    if (opt.Contains("--plotsdir=")){
      plotsDir = opt(11,opt.Sizeof());
    }
    if (opt.Contains("--maxResid=")){
      maxResid=((TString)opt(11,opt.Sizeof())).Atof();
    }
    if (opt.Contains("-l")){
      i++;
      label = argv[i];
    }
    if (opt == "--isMC"){
      isMC=1;
    }
  }
  


  gStyle->SetTitleSize(0.05,"XYZT");
  gStyle->SetPadLeftMargin(.15);
  gStyle->SetPadBottomMargin(.13);
  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);

  gStyle->SetPadLeftMargin(.10);
  gStyle->SetPadRightMargin(.075);
  TCanvas* c2 = new TCanvas("c2","c2",1200,900);
  
  gStyle->SetPadLeftMargin(.15);
  gStyle->SetPadRightMargin(.02);
  gStyle->SetPadBottomMargin(.12);
  TCanvas* c3 = new TCanvas("c3","c3",1200,1600);
  c3->Divide(3,4);
  
  TLegend* legend1 = new TLegend(0.15, 0.7, 0.90, 0.9);
  TLegend* legend2 = new TLegend(0.15, 0.7, 0.90, 0.9);
  TLegend* legend3 = new TLegend(0.15, 0.7, 0.90, 0.9);
  TLegend* legend4 = new TLegend(0.15, 0.7, 0.70, 0.9);
  
  TLegend* legend5 = new TLegend(0.1, 0.75, 0.40, 0.9);
  TLegend* legend6 = new TLegend(0.15, 0.70, 0.70, 0.9);
  
  
  
  for(int before_after = 0; before_after<2; before_after++){
    
    int n_chi2_lt_2 = 0;
    int n_tracks = 0;
    // auto save every MB
    //AlignTree->SetAutoSave(1000000);
    
    // these histograms are for debugging.  The limits and labels for
    // some of these histograms are taylored for the Clas12 CVT.
    
    /*if(before_after){
      gStyle->SetMarkerColor(kBlack);
      gStyle->SetLineColor(kBlack);
    } else {
      gStyle->SetMarkerColor(kRed);
      gStyle->SetLineColor(kRed);
    }*/
    
    TString suffix = before_after ? "After" : "Before";

    TH1F*  hchi2ndof = new TH1F("hchi2ndof"+suffix, "track #chi^{2}/ndof;track #chi^{2}/n_{dof};# of tracks", 100, 0, 20);
    
    TH1F* residuals_svt = new TH1F ("res_svt"+suffix, "SVT residuals;residual [mm];# of clusters", 100, -0.5, 0.5);
    TH1F* residuals_bmtz = new TH1F ("res_bmtz"+suffix, "BMTZ residuals;residual [mm];# of clusters", 100, -1.5, 1.5);
    TH1F* residuals_bmtc = new TH1F ("res_bmtc"+suffix, "BMTC residuals;residual [mm];# of clusters", 100, -1.5, 1.5);
    
    //shift the "before" by a fraction of a bin for clarity
    double shift_module = -(1.-before_after)*1./5;

    #define RESID_BINS 100, -3,3
    TH2D* residuals_vs_module =  new TH2D ("res_mod_"+suffix, "Residuals (all modules);module # ;residual [mm]", 102, 0+shift_module, 102+shift_module,
                                           RESID_BINS);
    
    //sets the error bar to be the std dev, instead of standard error on the mean
    //residuals_vs_module->SetErrorOption("s");
    
    int color = before_after ? kBlack: kRed;
     residuals_svt->SetLineColor(color);
     residuals_bmtz->SetLineColor(color);
     residuals_bmtc->SetLineColor(color);
     hchi2ndof->SetLineColor(color);
    
    residuals_vs_module->SetLineColor(color);
    residuals_vs_module->SetMarkerColor(color);
    
     int linestyle = before_after ? 1: 2;
     residuals_svt->SetLineStyle(linestyle);
     residuals_bmtz->SetLineStyle(linestyle);
     residuals_bmtc->SetLineStyle(linestyle);
     hchi2ndof->SetLineStyle(linestyle);
    
    int markerstyle=before_after ? 21: 25;
    residuals_vs_module->SetMarkerStyle(markerstyle);
    
    
    //number of bins for the residuals vs kinematics plots.
    int nbins_kin = 15;
    //shift the "before" markers by a 5th of a bin for visual clarity
    double shift_d  = -(1.-before_after)*(dmax-dmin)/nbins_kin/5;
    double shift_phi  = -(1.-before_after)*(phimax-phimin)/nbins_kin/5;
    double shift_z  = -(1.-before_after)*(zmax-zmin)/nbins_kin/5;
    double shift_tandip  = -(1.-before_after)*(tandipmax-tandipmin)/nbins_kin/5;
    cout << shift_d << endl;
    
    
    
    TH2D* residuals_vs_phi_svt = new  TH2D ("res_phi_svt"+suffix, "SVT residuals vs #phi_{0};#phi_{0} [rad];residual [mm]", nbins_kin,phimin, phimax, RESID_BINS);
    TH2D* residuals_vs_d0_svt =  new TH2D ("res_d0_svt"+suffix, "SVT residuals vs d_{0};d_{0} [mm]; residual [mm]", nbins_kin, dmin, dmax, RESID_BINS);
    TH2D* residuals_vs_theta_svt =  new TH2D ("res_theta_svt"+suffix, "SVT residuals vs t_{0};t_{0};residual [mm]", nbins_kin, tandipmin, tandipmax, RESID_BINS);
    TH2D* residuals_vs_z_svt =  new TH2D ("res_z_svt"+suffix, "SVT residuals vs z_{0};z_{0} [mm];residual [mm]", nbins_kin, zmin, zmax, RESID_BINS);
    
    TH2D* residuals_vs_phi_bmtz = new  TH2D ("res_phi_bmtz"+suffix, "BMTZ residuals vs #phi_{0};#phi_{0} [rad];residual [mm]", nbins_kin, phimin,phimax, RESID_BINS);
    TH2D* residuals_vs_d0_bmtz = new  TH2D ("res_d0_bmtz"+suffix, "BMTZ residuals vs d_{0};d_{0} [mm]; residual [mm]", nbins_kin, dmin, dmax, RESID_BINS);
    TH2D* residuals_vs_theta_bmtz = new  TH2D ("res_theta_bmtz"+suffix, "BMTZ residuals vs t_{0};t_{0};residual [mm]", nbins_kin, tandipmin, tandipmax, RESID_BINS);
    TH2D* residuals_vs_z_bmtz =  new TH2D ("res_z_bmtz"+suffix, "BMTZ residuals vs z_{0};z_{0} [mm];residual [mm]", nbins_kin, zmin, zmax, RESID_BINS);
    
    TH2D* residuals_vs_phi_bmtc = new  TH2D ("res_phi_bmtc"+suffix, "BMTC residuals vs #phi_{0};#phi_{0} [rad];residual [mm]", nbins_kin, phimin, phimax,RESID_BINS);
    TH2D* residuals_vs_d0_bmtc =  new TH2D ("res_d0_bmtc"+suffix, "BMTC residuals vs d_{0};d_{0} [mm]; residual [mm]", nbins_kin, dmin, dmax, RESID_BINS);
    TH2D* residuals_vs_theta_bmtc =  new TH2D ("res_theta_bmtc"+suffix, "BMTC residuals vs t_{0};t_{0};residual [mm]", nbins_kin, tandipmin, tandipmax, RESID_BINS);
    TH2D* residuals_vs_z_bmtc =  new TH2D ("res_z_bmtc"+suffix, "BMTC residuals vs z_{0};z_{0} [mm];residual [mm]", nbins_kin, zmin, zmax, RESID_BINS);
    
/*#define format_profile(profile,min,max) \
    
    
    format_profile(residuals_vs_d0_svt,-0.6,0.6)
    format_profile(residuals_vs_phi_svt,-0.6,0.6)
    format_profile(residuals_vs_z_svt,-0.6,0.6)
    format_profile(residuals_vs_theta_svt,-0.6,.6)
    format_profile(residuals_vs_d0_bmtz,-5,5)
    format_profile(residuals_vs_phi_bmtz,-5,5)
    format_profile(residuals_vs_z_bmtz,-5,5)
    format_profile(residuals_vs_theta_bmtz,-5,5)
    format_profile(residuals_vs_d0_bmtc,-5,5)
    format_profile(residuals_vs_phi_bmtc,-5,5)
    format_profile(residuals_vs_z_bmtc,-5,5)
    format_profile(residuals_vs_theta_bmtc,-5,5)
*/
    
    
    cout << "initialized histograms" << endl;
    int events = 0, tracks=0;
    
    int ndof =0;
    float chi2 =0;
    
    //if there is no input file
    
    TString inputFileName = before_after ? inputFileNameAfter : inputFileNameBefore;
    if(inputFileName==TString())  {
      std::cout << " *** please provide a file name..." << std::endl;
      exit(0);
    }
    
    
    
    /////////////////////////////////////
    
    TFile * inputFile = TFile::Open(inputFileName);
    cout << "found file"<<endl;
    
    AlignInfo * alignInfo = (AlignInfo *)inputFile->Get("AlignInfo");
    if(alignInfo == NULL){
      cout << "align info not found" << endl;
      exit(0);
    }
    cout << "align info"  << endl;
    int alignables = alignInfo->GetNAlignables();
    int parameters = alignInfo->GetNParameters();
    cout << alignables << "alignables; " <<  parameters << " parameters"<< endl;
    
    AlignEvent *aevent = new AlignEvent();
    TTree *AlignTree = (TTree*)inputFile->Get("AlignTree");
    AlignTree->SetBranchAddress("AlignEvent", & aevent);
    
    float d0;
    float phi;
    float z;
    float tandip;
    AlignTree->SetBranchAddress("q0", &d0);
    AlignTree->SetBranchAddress("q1", &phi);
    AlignTree->SetBranchAddress("q2", &z);
    AlignTree->SetBranchAddress("q3", &tandip);
    
    
    for(int i = 0; i<AlignTree->GetEntriesFast();i++){
      
      AlignTree->GetEntry(i);
      int ndof = aevent->GetNdof();
      float chi2 = aevent->GetChi2();
      
      //cout << d0 << " "<< phi << " "<< z << " " << tandip << endl;
      
      hchi2ndof->Fill(chi2/ndof);
      if(chi2/ndof<2)
        n_chi2_lt_2++;
      n_tracks++;
      
      
      for(int j = 0; j < aevent->GetMeasuredCovariance()->GetNrows(); j++){
        double resid = (*aevent->GetMeasurements())(j)-(*aevent->GetTrackPrediction())(j);
        double res = sqrt((*aevent->GetMeasuredCovariance())[j][j]);
        double module = aevent->GetIndex()->At(j);
        
        
        residuals_vs_module->Fill(module,resid);
       
        if(module >=84 && module < 102){
          if (module <=86 || module>=93 && module<=95 || module>=99){
            residuals_bmtc->Fill(resid);
            residuals_vs_d0_bmtc->Fill(d0, resid);
            residuals_vs_phi_bmtc->Fill(phi, resid);
            residuals_vs_z_bmtc->Fill(z, resid);
            residuals_vs_theta_bmtc->Fill(tandip, resid);
          } else {
            residuals_bmtz->Fill(resid);
            residuals_vs_d0_bmtz->Fill(d0, resid);
            residuals_vs_phi_bmtz->Fill(phi, resid);
            residuals_vs_z_bmtz->Fill(z, resid);
            residuals_vs_theta_bmtz->Fill(tandip, resid);
          }
        } else if (module == 102){
          //beamspot.
        } else {
          residuals_svt->Fill(resid);
          residuals_vs_d0_svt->Fill(d0, resid);
          residuals_vs_phi_svt->Fill(phi, resid);
          residuals_vs_z_svt->Fill(z, resid);
          residuals_vs_theta_svt->Fill(tandip, resid);
        }
        
      }
      
    }
    cout << suffix << endl;
    
    cout << suffix << ":  frac of events with chi2/ndof<2: "  << (n_chi2_lt_2/(double)n_tracks) << endl;
    
    TString opt = before_after ? "SAME" : "";
    
    if (isMC){
      residuals_svt->SetTitle(residuals_svt->GetTitle()+(TString)" (MC)");
      residuals_bmtz->SetTitle(residuals_bmtz->GetTitle()+(TString)" (MC)");
      residuals_bmtc->SetTitle(residuals_bmtc->GetTitle()+(TString)" (MC)");
      hchi2ndof->SetTitle(hchi2ndof->GetTitle()+(TString)" (MC)");
    }
    
    c1->cd(1);
    legend1->AddEntry(residuals_svt, Form(suffix + ",\n RMS = %.2f mm, fit #sigma = %.3f mm",residuals_svt->GetRMS(), getSigma(residuals_svt)),"l");
    residuals_svt->Draw(opt);
    residuals_svt->SetMaximum(residuals_svt->GetMaximum()*(isMC ? 7 : 4));
    c1->cd(2);
    legend2->AddEntry(residuals_bmtz, Form(suffix + ",\n RMS = %.2f mm, fit #sigma = %.2f mm",residuals_bmtz->GetRMS(), getSigma(residuals_bmtz)),"l");
    residuals_bmtz->Draw(opt);
    residuals_bmtz->SetMaximum(residuals_bmtz->GetMaximum()*(isMC ? 9 : 7));
    c1->cd(3);
    legend3->AddEntry(residuals_bmtc, Form(suffix+ ",\n RMS = %.2f mm, fit #sigma = %.2f mm",residuals_bmtc->GetRMS(), getSigma(residuals_bmtc)), "l");
    residuals_bmtc->Draw(opt);
    residuals_bmtc->SetMaximum(residuals_bmtc->GetMaximum()*(isMC ? 3 : 2.4));
    c1->cd(4);
    legend4->AddEntry(hchi2ndof, Form(suffix+", mean = %.1f",hchi2ndof->GetMean()),"l");
    hchi2ndof->SetMaximum(hchi2ndof->GetMaximum()*(isMC ? 21 : 10));
    hchi2ndof->Draw(opt);
    
   
    
    opt = before_after ? "SP" : "AP";
    
    TLine* line = new TLine();
    line->SetLineStyle(2);
    line->SetLineColorAlpha(kBlack,0.5);
    c2->cd();
    residuals_vs_module->GetYaxis()->SetTitleOffset(1.5);
    line->DrawLine(0,0,102,0);
    //residuals_vs_module->SetMinimum(-2);
    //residuals_vs_module->SetMaximum(3);
    if (isMC){
      suffix = "MC (" + suffix + ")";
      residuals_vs_module->SetTitle("Residuals (all modules, MC)");
    }
    TGraphErrors* t = createProfile(residuals_vs_module, color, markerstyle, shift_module,1.5);
    if(before_after){
      double worstMuSVT=0;
      double worstMuBMTZ=0;
      double worstMuBMTC=0;
      for (int k = 0; k<t->GetN(); k++){
        double mu = t->GetPointY(k);
        if(k >=84 && k < 102){
          if (k <=86 || k>=93 && k<=95 || k>=99){
            if (abs(mu)>worstMuBMTC)
              worstMuBMTC=abs(mu);
          } else {
            if (abs(mu)>worstMuBMTZ){
              worstMuBMTZ=abs(mu);
            }
          }
        } else {
          if(abs(mu)>worstMuSVT){
            worstMuSVT=abs(mu);
          }
        }
      }
      cout << "worst residuals for modules (SVT, BMTZ, BMTC): " << worstMuSVT <<", " << worstMuBMTZ << ", "<< worstMuBMTC<<  endl;
    }
    t->GetHistogram()->SetMinimum(-1.15);
    t->GetHistogram()->SetMaximum(1.7);
    t->GetHistogram()->GetYaxis()->SetTitleOffset(1);
    t->Draw(opt);
    
    legend5->AddEntry(residuals_vs_module,suffix, "lp");
    
    if (isMC){
      residuals_vs_d0_svt->SetTitle(residuals_vs_d0_svt->GetTitle()+(TString)" (MC)");
      residuals_vs_d0_bmtz->SetTitle(residuals_vs_d0_bmtz->GetTitle()+(TString)" (MC)");
      residuals_vs_d0_bmtc->SetTitle(residuals_vs_d0_bmtc->GetTitle()+(TString)" (MC)");
      
      residuals_vs_phi_svt->SetTitle(residuals_vs_phi_svt->GetTitle()+(TString)" (MC)");
      residuals_vs_phi_bmtz->SetTitle(residuals_vs_phi_bmtz->GetTitle()+(TString)" (MC)");
      residuals_vs_phi_bmtc->SetTitle(residuals_vs_phi_bmtc->GetTitle()+(TString)" (MC)");
      
      residuals_vs_z_svt->SetTitle(residuals_vs_z_svt->GetTitle()+(TString)" (MC)");
      residuals_vs_z_bmtz->SetTitle(residuals_vs_z_bmtz->GetTitle()+(TString)" (MC)");
      residuals_vs_z_bmtc->SetTitle(residuals_vs_z_bmtc->GetTitle()+(TString)" (MC)");
      
      residuals_vs_theta_svt->SetTitle(residuals_vs_theta_svt->GetTitle()+(TString)" (MC)");
      residuals_vs_theta_bmtz->SetTitle(residuals_vs_theta_bmtz->GetTitle()+(TString)" (MC)");
      residuals_vs_theta_bmtc->SetTitle(residuals_vs_theta_bmtc->GetTitle()+(TString)" (MC)");
    }
    
    
    c3->cd(1);TGraphErrors *graph = createProfile(residuals_vs_d0_svt, color, markerstyle, shift_d,0.5);graph->Draw(opt);
    legend6->AddEntry(graph,suffix, "lp");
    line->DrawLine(dmin,0,dmax,0);
    c3->cd(4);createProfile(residuals_vs_phi_svt, color, markerstyle, shift_phi,0.5)->Draw(opt);
    line->DrawLine(phimin,0,phimax,0);
    c3->cd(7);createProfile(residuals_vs_z_svt, color, markerstyle, shift_z,0.5)->Draw(opt);
    line->DrawLine(zmin,0,zmax,0);
    c3->cd(10);createProfile(residuals_vs_theta_svt, color, markerstyle, shift_tandip,0.5)->Draw(opt);
    line->DrawLine(tandipmin,0,tandipmax,0);
    
    
    c3->cd(2);createProfile(residuals_vs_d0_bmtz, color, markerstyle, shift_d)->Draw(opt);
    line->DrawLine(dmin,0,dmax,0);
    c3->cd(5);createProfile(residuals_vs_phi_bmtz, color, markerstyle, shift_phi)->Draw(opt);
    line->DrawLine(phimin,0,phimax,0);
    c3->cd(8);createProfile(residuals_vs_z_bmtz, color, markerstyle, shift_z)->Draw(opt);
    line->DrawLine(zmin,0,zmax,0);
    c3->cd(11);createProfile(residuals_vs_theta_bmtz, color, markerstyle, shift_tandip)->Draw(opt);
    line->DrawLine(tandipmin,0,tandipmax,0);
    
    c3->cd(3);createProfile(residuals_vs_d0_bmtc, color, markerstyle, shift_d)->Draw(opt);
    line->DrawLine(dmin,0,dmax,0);
    c3->cd(6);createProfile(residuals_vs_phi_bmtc, color, markerstyle, shift_phi)->Draw(opt);
    line->DrawLine(phimin,0,phimax,0);
    c3->cd(9);createProfile(residuals_vs_z_bmtc, color, markerstyle, shift_z)->Draw(opt);
    line->DrawLine(zmin,0,zmax,0);
    c3->cd(12);createProfile(residuals_vs_theta_bmtc, color, markerstyle, shift_tandip)->Draw(opt);
    //c3->cd(12);residuals_vs_theta_bmtc->Draw(opt);
    line->DrawLine(tandipmin,0,tandipmax,0);
    
    
    
    inputFile->Close();
  }
  
  
  TString ext = "pdf";
  if(!(plotsDir==TString())){
    c1->cd(1);
    legend1->Draw();
    c1->cd(2);
    legend2->Draw();
    c1->cd(3);
    legend3->Draw();
    c1->cd(4);
    legend4->Draw();
    c1->SaveAs(plotsDir+"/residuals_1d." +ext);
    
    c2->cd();
    legend5->Draw();
    TLine* line = new TLine();
    line->SetLineStyle(3);
    double xs[] = {42,84,87,93,96,99};
    for(double x : xs){
      if(x == 84)
        line->SetLineWidth(2);
      else
        line->SetLineWidth(1);
      line->DrawLine(x,-1.15,x,x <= 84? 1.7: 1.4);
    }
    
    TText *text = new TText();
    text->SetTextSize(.06);
    text->DrawText(12, 0.82, "SVT (inner)");
    text->DrawText(12+42, 0.82, "SVT (outer)");
    
    text->DrawText(90, 1.5, "BMT");
    text->SetTextSize(.04);
    
    text->DrawText(84+0.8, 1.3, "C");
    text->DrawText(87+0.8, 1.3, "Z");
    text->DrawText(90+0.8, 1.3, "Z");
    text->DrawText(93+0.8, 1.3, "C");
    text->DrawText(96+0.8, 1.3, "Z");
    text->DrawText(99+0.8, 1.3, "C");
    
    c2->SaveAs(plotsDir+"/residuals_module." + ext);
    c3->cd(1);
    legend6->Draw();
    //if(label != ""):
    //  text->DrawText(0, label);
    c3->SaveAs(plotsDir+"/residuals_kinematics." + ext);
  }
  
  
  //auto finish = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<double> elapsed = finish - start;
  //std::cout << "Elapsed time: " << elapsed.count()<< "s, events = "<<events<< "\n";
  
}
