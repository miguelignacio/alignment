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

using namespace std;

double dmin=-2,dmax = 2;
double phimin=-TMath::Pi(), phimax = TMath::Pi();
double tandipmin= -1,tandipmax =1.5;
double zmin=-80, zmax =10;

int main(int argc, char * argv[]) {
  // Record start time
  //auto start = std::chrono::high_resolution_clock::now();
  
  TString inputFileNameBefore;
  TString inputFileNameAfter;
  TString plotsDir;

  double maxResid=5;
  vector<TString> inputFiles;
  TString label ="";
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
  }
  


  gStyle->SetTitleSize(0.05,"XYZT");
  gStyle->SetPadLeftMargin(.15);
  gStyle->SetPadBottomMargin(.13);
  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);
  
  TCanvas* c2 = new TCanvas("c2","c2",1200,600);
  
  gStyle->SetPadLeftMargin(.15);
  gStyle->SetPadRightMargin(.02);
  gStyle->SetPadBottomMargin(.12);
  TCanvas* c3 = new TCanvas("c3","c3",1200,1600);
  c3->Divide(3,4);
  
  TLegend* legend1 = new TLegend(0.15, 0.7, 0.70, 0.9);
  TLegend* legend2 = new TLegend(0.15, 0.7, 0.70, 0.9);
  TLegend* legend3 = new TLegend(0.15, 0.7, 0.70, 0.9);
  TLegend* legend4 = new TLegend(0.15, 0.7, 0.70, 0.9);
  
  TLegend* legend5 = new TLegend(0.15, 0.7, 0.50, 0.9);
  
  
  
  for(int before_after = 0; before_after<2; before_after++){
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
    
    TH1F*  hchi2ndof = new TH1F("hchi2ndof"+suffix, "track #chi^{2}/ndof;track #chi^{2}/n_{dof};# of events", 100, 0, 40);
    
    TH1F* residuals_svt = new TH1F ("res_svt"+suffix, "SVT residuals;residual [mm];# of clusters", 100, -2, 2);
    TH1F* residuals_bmtz = new TH1F ("res_bmtz"+suffix, "BMTZ residuals;residual [mm];# of clusters", 100, -5, 5);
    TH1F* residuals_bmtc = new TH1F ("res_bmtc"+suffix, "BMTC residuals;residual [mm];# of clusters", 100, -5, 5);
    
   
    
    
    TProfile* residuals_vs_module =  new TProfile ("res_mod_"+suffix, "Residuals (all modules);module # ;residual [mm]", 103, 0, 103);
    
    //sets the error bar to be the std dev, instead of standard error on the mean
    residuals_vs_module->SetErrorOption("s");
    
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
    
    
    
    TProfile* residuals_vs_phi_svt = new  TProfile ("res_phi_svt"+suffix, "SVT residuals vs #phi;#phi [rad];residual [mm]", 30,phimin, phimax);
    TProfile* residuals_vs_d0_svt =  new TProfile ("res_d0_svt"+suffix, "SVT residuals vs d_{0};d_{0} [mm]; residual [mm]", 30, dmin, dmax);
    TProfile* residuals_vs_theta_svt =  new TProfile ("res_theta_svt"+suffix, "SVT residuals vs tan#theta_{dip};tan#theta_{dip};residual [mm]", 30, tandipmin, tandipmax);
    TProfile* residuals_vs_z_svt =  new TProfile ("res_z_svt"+suffix, "SVT residuals vs z;z [mm];residual [mm]", 30, zmin, zmax);
    
    TProfile* residuals_vs_phi_bmtz = new  TProfile ("res_phi_bmtz"+suffix, "BMTZ residuals vs #phi;#phi [rad];residual [mm]", 30,phimin,phimax);
    TProfile* residuals_vs_d0_bmtz = new  TProfile ("res_d0_bmtz"+suffix, "BMTZ residuals vs d_{0};d_{0} [mm]; residual [mm]", 30, dmin, dmax);
    TProfile* residuals_vs_theta_bmtz = new  TProfile ("res_theta_bmtz"+suffix, "BMTZ residuals vs tan#theta_{dip};tan#theta_{dip};residual [mm]", 30, tandipmin, tandipmax);
    TProfile* residuals_vs_z_bmtz =  new TProfile ("res_z_bmtz"+suffix, "BMTZ residuals vs z;z [mm];residual [mm]", 30, zmin, zmax);
    
    TProfile* residuals_vs_phi_bmtc = new  TProfile ("res_phi_bmtc"+suffix, "BMTC residuals vs #phi;#phi [rad];residual [mm]", 30, phimin, phimax);
    TProfile* residuals_vs_d0_bmtc =  new TProfile ("res_d0_bmtc"+suffix, "BMTC residuals vs d_{0};d_{0} [mm]; residual [mm]", 30, dmin, dmax);
    TProfile* residuals_vs_theta_bmtc =  new TProfile ("res_theta_bmtc"+suffix, "BMTC residuals vs tan#theta_{dip};tan#theta_{dip};residual [mm]", 30, tandipmin, tandipmax);
    TProfile* residuals_vs_z_bmtc =  new TProfile ("res_z_bmtc"+suffix, "BMTC residuals vs z;z [mm];residual [mm]", 30, zmin, zmax);
    
#define format_profile(profile,min,max) \
    profile->SetMarkerStyle(markerstyle);\
    profile->SetMarkerColor(color);\
    profile->SetLineColor(color);\
    profile->SetMinimum(min);\
    profile->SetMaximum(max);\
    profile->SetErrorOption("s");
    
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
    
    /*
    residuals_vs_d0_svt->SetMarkerStyle(markerstyle);
    residuals_vs_phi_svt->SetMarkerStyle(markerstyle);
    residuals_vs_theta_svt->SetMarkerStyle(markerstyle);
    residuals_vs_z_svt->SetMarkerStyle(markerstyle);
    residuals_vs_d0_bmtz->SetMarkerStyle(markerstyle);
    residuals_vs_phi_bmtz->SetMarkerStyle(markerstyle);
    residuals_vs_theta_bmtz->SetMarkerStyle(markerstyle);
    residuals_vs_z_bmtz->SetMarkerStyle(markerstyle);
    residuals_vs_d0_bmtc->SetMarkerStyle(markerstyle);
    residuals_vs_phi_bmtc->SetMarkerStyle(markerstyle);
    residuals_vs_theta_bmtc->SetMarkerStyle(markerstyle);
    residuals_vs_z_bmtc->SetMarkerStyle(markerstyle);
    
    residuals_vs_d0_svt->SetMarkerColor(color);
    residuals_vs_phi_svt->SetMarkerColor(color);
    residuals_vs_theta_svt->SetMarkerColor(color);
    residuals_vs_z_svt->SetMarkerColor(color);
    residuals_vs_d0_bmtz->SetMarkerColor(color);
    residuals_vs_phi_bmtz->SetMarkerColor(color);
    residuals_vs_theta_bmtz->SetMarkerColor(color);
    residuals_vs_z_bmtz->SetMarkerColor(color);
    residuals_vs_d0_bmtc->SetMarkerColor(color);
    residuals_vs_phi_bmtc->SetMarkerColor(color);
    residuals_vs_theta_bmtc->SetMarkerColor(color);
    residuals_vs_z_bmtc->SetMarkerColor(color);
    
    residuals_vs_d0_svt->SetLineColor(color);
    residuals_vs_phi_svt->SetLineColor(color);
    residuals_vs_theta_svt->SetLineColor(color);
    residuals_vs_z_svt->SetLineColor(color);
    residuals_vs_d0_bmtz->SetLineColor(color);
    residuals_vs_phi_bmtz->SetLineColor(color);
    residuals_vs_theta_bmtz->SetLineColor(color);
    residuals_vs_z_bmtz->SetLineColor(color);
    residuals_vs_d0_bmtc->SetLineColor(color);
    residuals_vs_phi_bmtc->SetLineColor(color);
    residuals_vs_theta_bmtc->SetLineColor(color);
    residuals_vs_z_bmtc->SetLineColor(color);
    
    residuals_vs_d0_svt->SetMaximum(0.3);
    residuals_vs_phi_svt->SetMaximum(0.3);
    residuals_vs_theta_svt->SetMaximum(0.3);
    residuals_vs_z_svt->SetMaximum(0.3);
    residuals_vs_d0_bmtz->SetMaximum(3);
    residuals_vs_phi_bmtz->SetMaximum(3);
    residuals_vs_theta_bmtz->SetMaximum(3);
    residuals_vs_z_bmtz->SetMaximum(3);
    residuals_vs_d0_bmtc->SetMaximum(3);
    residuals_vs_phi_bmtc->SetMaximum(3);
    residuals_vs_theta_bmtc->SetMaximum(3);
    residuals_vs_z_bmtc->SetMaximum(3);
    
    residuals_vs_d0_svt->SetMinimum(-0.3);
    residuals_vs_phi_svt->SetMinimum(-0.3);
    residuals_vs_theta_svt->SetMinimum(-0.3);
    residuals_vs_z_svt->SetMinimum(-0.3);
    residuals_vs_d0_bmtz->SetMinimum(-3);
    residuals_vs_phi_bmtz->SetMinimum(-3);
    residuals_vs_theta_bmtz->SetMinimum(-3);
    residuals_vs_z_bmtz->SetMinimum(-3);
    residuals_vs_d0_bmtc->SetMinimum(-3);
    residuals_vs_phi_bmtc->SetMinimum(-3);
    residuals_vs_theta_bmtc->SetMinimum(-3);
    residuals_vs_z_bmtc->SetMinimum(-3);*/
    
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
      
      cout << d0 << " "<< phi << " "<< z << " " << tandip << endl;
      
      hchi2ndof->Fill(chi2/ndof);
      
      
      
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
    
    
    
    TString opt = before_after ? "SAME" : "";
    c1->cd(1);
    legend1->AddEntry(residuals_svt, Form(suffix + ",\n RMS = %.1f mm",residuals_svt->GetRMS()),"l");
    residuals_svt->Draw(opt);
    residuals_svt->SetMaximum(residuals_svt->GetMaximum()*2);
    c1->cd(2);
    legend2->AddEntry(residuals_bmtz, Form(suffix + ",\n RMS = %.1f mm",residuals_bmtz->GetRMS()),"l");
    residuals_bmtz->Draw(opt);
    residuals_bmtz->SetMaximum(residuals_bmtz->GetMaximum()*3);
    c1->cd(3);
    legend3->AddEntry(residuals_bmtc, Form(suffix+ ",\n RMS = %.1f mm",residuals_bmtc->GetRMS()), "l");
    residuals_bmtc->Draw(opt);
    residuals_bmtc->SetMaximum(residuals_bmtc->GetMaximum()*2);
    c1->cd(4);
    legend4->AddEntry(hchi2ndof, Form(suffix+", mean = %.1f",hchi2ndof->GetMean()),"l");
    hchi2ndof->SetMaximum(hchi2ndof->GetMaximum()*3);
    hchi2ndof->Draw(opt);
    
    
    c2->cd();
    residuals_vs_module->Draw(opt);
    legend5->AddEntry(residuals_vs_module,suffix, "lp");
    
    
    TLine* line = new TLine();
    line->SetLineStyle(2);
    c3->cd(1);residuals_vs_d0_svt->Draw(opt);
    line->DrawLine(dmin,0,dmax,0);
    c3->cd(4);residuals_vs_phi_svt->Draw(opt);
    line->DrawLine(phimin,0,phimax,0);
    c3->cd(7);residuals_vs_z_svt->Draw(opt);
    line->DrawLine(zmin,0,zmax,0);
    c3->cd(10);residuals_vs_theta_svt->Draw(opt);
    line->DrawLine(tandipmin,0,tandipmax,0);
    
    
    c3->cd(2);residuals_vs_d0_bmtz->Draw(opt);
    line->DrawLine(dmin,0,dmax,0);
    c3->cd(5);residuals_vs_phi_bmtz->Draw(opt);
    line->DrawLine(phimin,0,phimax,0);
    c3->cd(8);residuals_vs_z_bmtz->Draw(opt);
    line->DrawLine(zmin,0,zmax,0);
    c3->cd(11);residuals_vs_theta_bmtz->Draw(opt);
    line->DrawLine(zmin,0,zmax,0);
    
    c3->cd(3);residuals_vs_d0_bmtc->Draw(opt);
    line->DrawLine(dmin,0,dmax,0);
    c3->cd(6);residuals_vs_phi_bmtc->Draw(opt);
    line->DrawLine(phimin,0,phimax,0);
    c3->cd(9);residuals_vs_z_bmtc->Draw(opt);
    line->DrawLine(zmin,0,zmax,0);
    c3->cd(12);residuals_vs_theta_bmtc->Draw(opt);
    line->DrawLine(tandipmin,0,tandipmax,0);
    
    
    
    inputFile->Close();
  }
  
  
  TString ext = "png";
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
    c2->SaveAs(plotsDir+"/residuals_module." + ext);
    c3->cd(1);
    legend5->Draw();
    c3->SaveAs(plotsDir+"/residuals_kinematics." + ext);
  }
  
  
  //auto finish = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<double> elapsed = finish - start;
  //std::cout << "Elapsed time: " << elapsed.count()<< "s, events = "<<events<< "\n";
  
}
