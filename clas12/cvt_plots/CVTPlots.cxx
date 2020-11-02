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
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TStyle.h>
#include <TLine.h>
#include <TText.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <TLegend.h>
#include "../event/AlignEvent.h"


#include <iostream>



#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"

using namespace std;


TH2F* createPlotXY(TH2F* h, TString title, TString ztitle){
  
  double R = 150;
  int Nbins = 200;
  TH2F*  hxy =   new TH2F(h->GetName()+(TString)"_xy", title +";x (mm);y (mm);"+ztitle, Nbins,-R,R,Nbins,-R,R);

  //cout << prof->GetNbinsX() << "bins" << endl;
  for(int mm = 0; mm<92; mm++){
    int m = mm%50;
    double n = h->ProjectionX("temp",mm+1,mm+2)->GetMean();
    
    double pi = TMath::Pi();
    double phim =  m<10 ? 2*pi*m/10 : (m<24 ? 2*pi*(m-10)/14 : 2*pi*(m-24)/18);
    double dphim = m<10 ? pi/10 : (m<24 ? pi/14 : pi/18);
    double rm = m<10 ? 65 : (m<24 ? 93 : 120);
    phim -= pi/2;
    while(phim < -pi)
      phim += 2*pi;
    while(phim > pi)
      phim -= 2*pi;
    for(int i = 0; i<hxy->GetNbinsX(); i++){
      double x = -R +2*i*R/Nbins;
      for(int j = 0; j<hxy->GetNbinsY(); j++){
        double y = -R +2*j*R/Nbins;
        double phi = TMath::ATan2(y,x);
        double rc = x*cos(phim)+y*sin(phim);
        if(mm<42 && rc>rm && rc<rm+2.3 && abs(phi-phim)<dphim){
          hxy->SetBinContent(i,j,n);
        }
	else if(mm>=50 && rc>rm+2.7 && rc<rm+5.0 && abs(phi-phim)<dphim){
          hxy->SetBinContent(i,j,n);
        }   
      }
    }
  }
  hxy->SetStats(0);
  
  return hxy;
}


int main(int argc, char * argv[]) {
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();
  
  TString inputFileName;
  TString plotsDir;
  
  vector<TString> inputFiles;
  for(Int_t i=1;i<argc;i++){
    TString opt=argv[i];
    if((opt.Contains("--in="))){
      inputFileName=opt(5,opt.Sizeof());

    } 
      
    if (opt.Contains("--plotsdir=")){
      plotsDir = opt(11,opt.Sizeof());
    }
  }
  //if there is no input file
  if(inputFileName==TString())  {
    std::cout << " *** please provide a file name..." << std::endl;
    exit(0);
  }
  /////////////////////////////////////
  
  TFile * inputFile = TFile::Open(inputFileName);
  
    
  cout << "found file"<<endl;


  AlignInfo * alignInfo = (AlignInfo *)inputFile->Get("AlignInfo");
  cout << "align info"  << endl;
  int alignables = alignInfo->GetNAlignables();
  int parameters = alignInfo->GetNParameters();
  cout << alignables << "alignables; " <<  parameters << " parameters"<< endl;
  
  AlignEvent *aevent = new AlignEvent();
  TTree *AlignTree = (TTree*)inputFile->Get("AlignTree");
  AlignTree->SetBranchAddress("AlignEvent", & aevent);
  // auto save every MB
  //AlignTree->SetAutoSave(1000000);

  // these histograms are for debugging.  The limits and labels for
  // some of these histograms are taylored for the Clas12 CVT.
  gStyle->SetTitleSize(0.05,"XYZT");
  gStyle->SetPadLeftMargin(.13);
  gStyle->SetPadBottomMargin(.13);

  TH1F*  hchi2 = new TH1F("hchi2", "track #chi^{2};track #chi^{2};# of events", 100, 0, 300);
  TH1F*  hchi2ndof = new TH1F("hchi2ndof", "track #chi^{2}/ndof;track #chi^{2}/n_{dof};# of events", 100, 0, 70);
  TH1F*  hchi2prob = new TH1F("hchi2prob", "track #chi^{2} probability;p(#chi^2,n_{dof};# of events", 100, 0, 1);
  TH1F*  hndof = new TH1F("hndof", "degrees of freedom;n_{dof};# of events", 20, 0, 20);
  TH1F*  hncross = new TH1F("hncross", "number of crosses used;# of crosses;# of events", 20, 0, 20);
  TH1F*  hres = new TH1F("hres", "resolution;resolution (mm);# of clusters", 100, 0, 0.1);
  TH2F*  hresmodule =	new TH2F("hresmodule", "resolution (mm);resolution;module id", 100, 0, 0.1,100,0,100);
  TH1F*  hmodule = new TH1F("hmodule", "module;module id;# of crosses", 50, 0, alignables ==42? 50: 100);
  
  TH1F*  htrackderiv1 = new TH1F("htrackderiv1", "Track derivatives (doca);B_{i,0} (dimensionless);# of clusters", 100, -2, 2);
  TH1F*  htrackderiv2 = new TH1F("htrackderiv2", "Track derivatives (phi0);B_{i,1} (mm);# of clusters", 100, -300, 300);
  TH1F*  htrackderiv3 = new TH1F("htrackderiv3", "Track derivatives (z0);B_{i,2} (dimensionless);# of clusters", 100, -0.1, 0.1);
  TH1F*  htrackderiv4 = new TH1F("htrackderiv4", "Track derivatives (tan dep);B_{i,3} (mm);# of clusters", 100, -30, 30);
  TH2F*  htrackderiv1mod = new TH2F("htrackderiv1mod", "Track derivatives (doca);B_{i,0} (dimensionless);module id", 100, -4, 4,100,0,100);
  TH2F*  htrackderiv2mod = new TH2F("htrackderiv2mod", "Track derivatives (phi0);B_{i,1} (mm);module id", 100, -300, 300,100,0,100);
  TH2F*  htrackderiv3mod = new TH2F("htrackderiv3mod", "Track derivatives (z0);B_{i,2} (dimensionless);module id", 100, -0.1, 0.1,100,0,100);
  TH2F*  htrackderiv4mod = new TH2F("htrackderiv4mod", "Track derivatives (tan dep);B_{i,3} (mm);module id", 100, -30, 30,100,0,100);

  TH1F*  halignderiv1 = new TH1F("halignderiv1", "Alignment derivatives (shift x);A_{i,0} (dimensionless);# of clusters", 100, -4, 4);
  TH1F*  halignderiv2 = new TH1F("halignderiv2", "Alignment derivatives (shift y);A_{i,1} (dimensionless);# of clusters", 100, -4, 4);
  TH1F*  halignderiv3 = new TH1F("halignderiv3", "Alignment derivatives (shift z);A_{i,2} (dimensionless);# of clusters", 100, -0.03, 0.07);
  TH1F*  halignderiv4 = new TH1F("halignderiv4", "Alignment derivatives (rotate x);A_{i,3} (mm/rad);# of clusters", 100, -300, 300);
  TH1F*  halignderiv5 = new TH1F("halignderiv5", "Alignment derivatives (rotate y);A_{i,4} (mm/rad);# of clusters", 100, -300, 300);
  TH1F*  halignderiv6 = new TH1F("halignderiv6", "Alignment derivatives (rotate z);A_{i,5} (mm/rad);# of clusters", 100, -300, 300);

  TH2F*  halignderiv1mod = new TH2F("halignderiv1mod", "Alignment derivatives (shift x);A_{i,0} (dimensionless);module id", 100, \
-4, 4,100,0,100);
  TH2F*  halignderiv2mod = new TH2F("halignderiv2mod", "Alignment derivatives (shift y);A_{i,1} (dimensionless);module id", 100, \
-1, 1, 100,0,100);
  TH2F*  halignderiv3mod = new TH2F("halignderiv3mod", "Alignment derivatives (shift z);A_{i,2} (dimensionless);module id", 100, \
-0.03, 0.07,100,0,100);
  TH2F*  halignderiv4mod = new TH2F("halignderiv4mod", "Alignment derivatives (rotate x);A_{i,3} (mm/rad);module id", 100, \
-300, 300,100,0,100);
  TH2F*  halignderiv5mod = new TH2F("halignderiv5mod", "Alignment derivatives (rotate y);A_{i,4} (mm/rad);module id", 100, \
-300, 300,100,0,100);
  TH2F*  halignderiv6mod = new TH2F("halignderiv6mod", "Alignment derivatives (rotate z);A_{i,5} (mm/rad);module id", 100, \
-300, 300,100,0,100);

  TH1F*  hmeas = new TH1F("hmeas", "1D position (measured);1D position (mm);# of clusters", 100, -50, 50);
  TH1F*  hextrap = new TH1F("hextrap", "1D position (extrap);1D position (mm);# of clusters", 100, -50, 50);
  TH2F*  hmeasextrap = new TH2F("hmeasextrap", "1D position (extrap);meas. 1D position (mm);extrap. 1D position (mm)", 100, -50, 50,100,-50,50);

  TH1F*  hmeasres = new TH1F("hmeasres", "1D position residual;meas - extrap (mm);# of clusters", 100, -2, 2);
  TH2F*  hmeasresmod = new TH2F("hmeasresmod", "1D position residual;meas - extrap (mm);module #", 100, -2, 2,100,0,100);

  cout << "initialized histograms" << endl;
  int events = 0, tracks=0;

  int ndof =0;
  float chi2 =0;
  for(int i = 0; i<AlignTree->GetEntriesFast();i++){

    AlignTree->GetEntry(i);
    int ndof = aevent->GetNdof();
    float chi2 = aevent->GetChi2();

    
    hchi2->Fill(chi2);


    hchi2ndof->Fill(chi2/ndof);
    hchi2prob->Fill(TMath::Prob(chi2, ndof));
    hndof->Fill(ndof);
    hncross->Fill(aevent->GetMeasuredCovariance()->GetNrows()/2);
    

    TMatrixD & A = *aevent->GetAlignmentDerivatives();
    /*
    cout << "A" << endl;
    A.Print("");
    cout << "index" << endl;
    for (int j = 0; j < aevent->GetMeasuredCovariance()->GetNrows(); j++){
      cout <<  aevent->GetIndex()->At(j) << "\t" ;
    }
    cout << endl<< endl;
    */
    TMatrixD & B = *aevent->GetTrackDerivatives();

    for(int j = 0; j < aevent->GetMeasuredCovariance()->GetNrows(); j++){
      double res = sqrt((*aevent->GetMeasuredCovariance())[j][j]);
      double module = aevent->GetIndex()->At(j);
      hres->Fill(res);
      if(j%2 == 0 || alignables == 84) //only one per cross
	hmodule->Fill(module);
      if(alignables == 42)
	module+= 50*(j%2);
      else
	module+= 8*(j%2); //the xy plotter assumes that the outer modules indices are offeset by 50 instead of 42
      hresmodule->Fill(res,module);


      

      
      htrackderiv1->Fill(B[j][0]);
      htrackderiv2->Fill(B[j][1]);
      htrackderiv3->Fill(B[j][2]);
      htrackderiv4->Fill(B[j][3]);
      htrackderiv1mod->Fill(B[j][0],module);
      htrackderiv2mod->Fill(B[j][1],module);
      htrackderiv3mod->Fill(B[j][2],module);
      htrackderiv4mod->Fill(B[j][3],module);

      int offset = alignables == 42 ? (j/2)*parameters : j*parameters;
      halignderiv1->Fill(A[j][offset+0]);
      if(parameters>1)
      halignderiv2->Fill(A[j][offset+1]);
      if(parameters>2)
      halignderiv3->Fill(A[j][offset+2]);
      if(parameters>3)
      halignderiv4->Fill(A[j][offset+3]);
      if(parameters>4)
      halignderiv5->Fill(A[j][offset+4]);
      if(parameters>5)
      halignderiv6->Fill(A[j][offset+5]);
      
      halignderiv1mod->Fill(A[j][offset+0],module);
      if(parameters>1)
      halignderiv2mod->Fill(A[j][offset+1],module);
      if(parameters>2)
      halignderiv3mod->Fill(A[j][offset+2],module);
      if(parameters>3)
      halignderiv4mod->Fill(A[j][offset+3],module);
      if(parameters>4)
      halignderiv5mod->Fill(A[j][offset+4],module);
      if(parameters>5)
      halignderiv6mod->Fill(A[j][offset+5],module);
      hmeas->Fill((*aevent->GetMeasurements())(j));
      hextrap->Fill((*aevent->GetTrackPrediction())(j));
      hmeasextrap->Fill((*aevent->GetMeasurements())(j),(*aevent->GetTrackPrediction())(j));
      hmeasres->Fill((*aevent->GetMeasurements())(j)-(*aevent->GetTrackPrediction())(j));
      hmeasresmod->Fill((*aevent->GetMeasurements())(j)-(*aevent->GetTrackPrediction())(j),module);
    }
  }
  
  //this is specific to the CVT
  double R = 150;
  int Nbins = 200;
  TH2F*  hmodulexy =   new TH2F("hresmodulexy", "modules xy;x (mm);y (mm)", Nbins,-R,R,Nbins,-R,R);
  for(int m = 0; m<42; m++){
    
    double n = hmodule->GetBinContent(m+1);
    double pi = TMath::Pi();
    double phim =  m<10 ? 2*pi*m/10 : (m<24 ? 2*pi*(m-10)/14 : 2*pi*(m-24)/18);
    double dphim = m<10 ? pi/10 : (m<24 ? pi/14 : pi/18);
    double rm = m<10 ? 65 : (m<24 ? 93 : 120);
    phim -= pi/2;
    while(phim < -pi)
      phim += 2*pi;
    while(phim > pi)
      phim -= 2*pi;
    for(int i = 0; i<hmodulexy->GetNbinsX(); i++){
      double x = -R +2*i*R/Nbins;
      for(int j = 0; j<hmodulexy->GetNbinsY(); j++){
	double y = -R +2*j*R/Nbins;
	double phi = TMath::ATan2(y,x);
	double rc = x*cos(phim)+y*sin(phim);
	if(rc>rm && rc<rm+5 && abs(phi-phim)<dphim){
	  hmodulexy->SetBinContent(i,j,n);
	}
      }
    }
  }
  
  
  if(!(plotsDir==TString())){
    TLine* l= new TLine(); l->SetLineColor(kRed);
    TCanvas* c = new TCanvas("canvas","canvas",800,600);
    hchi2->Draw();c->SaveAs(plotsDir+"/chi2.pdf");
    hchi2ndof->Draw();c->SaveAs(plotsDir+"/chi2ndof.pdf");
    hchi2prob->Draw();c->SaveAs(plotsDir+"/chi2prob.pdf");
    hndof->Draw();c->SaveAs(plotsDir+"/ndof.pdf");
    hncross->Draw();c->SaveAs(plotsDir+"/ncross.pdf");
    hres->Draw();c->SaveAs(plotsDir+"/resolution.pdf");
    hresmodule->Draw("COLZ1");c->SaveAs(plotsDir+"/res_vs_module.pdf");
    
    hmodule->Draw();

    l->DrawLine(10,0,10,hmodule->GetMaximum());
    l->DrawLine(24,0,24,hmodule->GetMaximum());
    l->DrawLine(42,0,42,hmodule->GetMaximum());
    TText * t = new TText();
    t->SetTextColor(kRed);
    t->DrawText(3,hmodule->GetMaximum()/10,"R1");
    t->DrawText(15,hmodule->GetMaximum()/10,"R2");
    t->DrawText(31,hmodule->GetMaximum()/10,"R3");
    c->SaveAs(plotsDir+"/module.pdf");
    
    htrackderiv1->Draw();c->SaveAs(plotsDir+"/B1.pdf");
    htrackderiv2->Draw();c->SaveAs(plotsDir+"/B2.pdf");
    htrackderiv3->Draw();c->SaveAs(plotsDir+"/B3.pdf");
    htrackderiv4->Draw();c->SaveAs(plotsDir+"/B4.pdf");
    htrackderiv1mod->Draw("COLZ1");c->SaveAs(plotsDir+"/B1mod.pdf");
    htrackderiv2mod->Draw("COLZ1");c->SaveAs(plotsDir+"/B2mod.pdf");
    htrackderiv3mod->Draw("COLZ1");c->SaveAs(plotsDir+"/B3mod.pdf");
    htrackderiv4mod->Draw("COLZ1");c->SaveAs(plotsDir+"/B4mod.pdf");
    gStyle->SetPadRightMargin(.2);
    createPlotXY(htrackderiv1mod,"B1","B1 (dimensionless)")->Draw("COLZ1");c->SaveAs(plotsDir+"/B1xy.pdf");
    createPlotXY(htrackderiv2mod,"B2","B2 (dimensionless)")->Draw("COLZ1");c->SaveAs(plotsDir+"/B2xy.pdf");
    createPlotXY(htrackderiv3mod,"B3","B3 (mm)")->Draw("COLZ1");c->SaveAs(plotsDir+"/B3xy.pdf");
    createPlotXY(htrackderiv4mod,"B4","B4 (mm)")->Draw("COLZ1");c->SaveAs(plotsDir+"/B4xy.pdf");
    gStyle->SetPadRightMargin(.10);
    halignderiv1->Draw();c->SaveAs(plotsDir+"/A1.pdf");
    halignderiv2->Draw();c->SaveAs(plotsDir+"/A2.pdf");
    halignderiv3->Draw();c->SaveAs(plotsDir+"/A3.pdf");
    halignderiv4->Draw();c->SaveAs(plotsDir+"/A4.pdf");
    halignderiv5->Draw();c->SaveAs(plotsDir+"/A5.pdf");
    halignderiv6->Draw();c->SaveAs(plotsDir+"/A6.pdf");
    halignderiv1mod->Draw("COLZ1");c->SaveAs(plotsDir+"/A1mod.pdf");
    halignderiv2mod->Draw("COLZ1");c->SaveAs(plotsDir+"/A2mod.pdf");
    halignderiv3mod->Draw("COLZ1");c->SaveAs(plotsDir+"/A3mod.pdf");
    halignderiv4mod->Draw("COLZ1");c->SaveAs(plotsDir+"/A4mod.pdf");
    halignderiv5mod->Draw("COLZ1");c->SaveAs(plotsDir+"/A5mod.pdf");
    halignderiv6mod->Draw("COLZ1");c->SaveAs(plotsDir+"/A6mod.pdf");

    gStyle->SetPadRightMargin(.2);
    createPlotXY(halignderiv1mod,"A1","A1 (dimensionless)")->Draw("COLZ1"); c->SaveAs(plotsDir+"/A1xy.pdf");
    createPlotXY(halignderiv2mod,"A2","A2 (dimensionless)")->Draw("COLZ1"); c->SaveAs(plotsDir+"/A2xy.pdf");
    createPlotXY(halignderiv3mod,"A3","A3 (dimensionless)")->Draw("COLZ1"); c->SaveAs(plotsDir+"/A3xy.pdf");
    createPlotXY(halignderiv4mod,"A4","A4 (mm)")->Draw("COLZ1"); c->SaveAs(plotsDir+"/A4xy.pdf");
    createPlotXY(halignderiv5mod,"A5","A5 (mm)")->Draw("COLZ1"); c->SaveAs(plotsDir+"/A5xy.pdf");
    createPlotXY(halignderiv6mod,"A6","A6 (mm)")->Draw("COLZ1"); c->SaveAs(plotsDir+"/A6xy.pdf");
    gStyle->SetPadRightMargin(.10);
    
    hmodulexy->Draw("COLZ1");c->SaveAs(plotsDir+"/modulexy.pdf");
    
    hmeas->SetTitle("1D hit positions");
    hmeas->Draw();hextrap->SetLineColor(kRed);hextrap->Draw("SAME");
    TLegend *legend = new TLegend();
    legend->AddEntry(hmeas, "measured");
    legend->AddEntry(hextrap, "track extrap");
    legend->Draw();
    c->SaveAs(plotsDir+"/meas_extrap.pdf");
    
    hmeasextrap->Draw("COLZ1");l->DrawLine(-50,-50,50,50);c->SaveAs(plotsDir+"/meas_vs_extrap.pdf");
    hmeasres->Draw();l->DrawLine(0,0,0,hmeasres->GetMaximum());c->SaveAs(plotsDir+"/measres.pdf");
    hmeasresmod->Draw("COLZ1");l->DrawLine(0,0,0,100);c->SaveAs(plotsDir+"/measresmod.pdf");

    TH1F*  hmeasres_shift = new TH1F("hmeasres_shift", "1D position residual;meas - extrap (mm);# of clusters", 100, -2, 2);
  
    
    gStyle->SetPadRightMargin(.2);
    createPlotXY(hmeasresmod,"residuals","residual (mm)")->Draw("COLZ1"); c->SaveAs(plotsDir+"/measresxy.pdf");
    
  }
  
  
   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< "s, events = "<<events<< "\n";
  
}
