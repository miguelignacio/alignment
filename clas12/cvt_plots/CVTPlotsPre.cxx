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

double angle(double phi){
  double pi = TMath::Pi();
  while(phi>pi)
    phi-=2*pi;
  while(phi<-pi)
    phi+=2*pi;
  return phi;
}

void drawLabel(TString& label){
  //if(label == "")
  //  label = "unlabeled";
  TText *pt = new TText(.15,.7, label);
  pt->SetNDC();
  pt->Draw();
}

double getMeanProjectionX(TH2* h, int j){
  double sumX = 0, sum = 0;
  for (int i = 0; i<h->GetNbinsX(); i++){
    double x = h->GetXaxis()->GetBinCenter(i+1);
    double bc = h->GetBinContent(i+1,j+1);
    sumX+= x*bc;
    sum+= bc;
  }
  return sumX/sum;
}
//change mode to "count" to make it show the number of hits on a given module
TH2F* createPlotXY(TH1* h, TString title, TString ztitle, TString mode ="mean"){
  
  double R = 235;
  int Nbins = 235;
  TH2F*  hxy =   new TH2F(h->GetName()+(TString)"_xy", title +";x (mm);y (mm);"+ztitle, Nbins,-R,R,Nbins,-R,R);

  //cout << prof->GetNbinsX() << "bins" << endl;
  for(int mm = 0; mm<102; mm++){
    double n = mode.EqualTo("mean") ? getMeanProjectionX((TH2*)h,mm) : h->GetBinContent(mm+1);
    if(mm < 84){
      int m = mm%42;
      if(TMath::IsNaN(n))
	continue;

      
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
	  else if(mm>=42 && rc>rm+2.7 && rc<rm+5.0 && abs(phi-phim)<dphim){
	    hxy->SetBinContent(i,j,n);
	  }   
	}
      }
    } else if(mm>=84){ //BMT:
      int m = mm-84;

      int layer = m/3;
      int sector = m%3;

      double radii[] = {146,161,176,191,206,221};
      double rm = radii[layer], drm = 4;
      double phisector = -TMath::Pi()/6+sector*TMath::Pi()*2/3;
      //cout<<  mm << " " << layer << " " << sector << " " << n << " " << phisector*180/TMath::Pi() << endl;
      
      for(int i = 0; i<hxy->GetNbinsX(); i++){
	double x = -R +2*i*R/Nbins;
	for(int j = 0; j<hxy->GetNbinsY(); j++){
	  double y = -R +2*j*R/Nbins;
	  double r = TMath::Hypot(x,y);
	  double phi = TMath::ATan2(y,x);
	  double dphi = TMath::Pi()/3-0.1;//the width between sectors is not right, but whatever.  
	  if(r>rm && r<=rm+drm && abs(angle(phi-phisector))<dphi){
	    hxy->SetBinContent(i,j,n);
	  }
	}
      }
    }
  }
  h->SetStats(0);
  hxy->SetStats(0);
  
  return hxy;
}


int main(int argc, char * argv[]) {
  // Record start time
  //auto start = std::chrono::high_resolution_clock::now();
  
  TString inputFileName;
  TString plotsDir;

  double maxResid=5;
  vector<TString> inputFiles;
  TString label ="";
  for(Int_t i=1;i<argc;i++){
    TString opt=argv[i];
    if((opt.Contains("--in="))){
      inputFileName=opt(5,opt.Sizeof());

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
  //if there is no input file
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
  // auto save every MB
  //AlignTree->SetAutoSave(1000000);

  // these histograms are for debugging.  The limits and labels for
  // some of these histograms are taylored for the Clas12 CVT.
  gStyle->SetTitleSize(0.05,"XYZT");
  gStyle->SetPadLeftMargin(.13);
  gStyle->SetPadBottomMargin(.13);

  TH1F*  hchi2 = new TH1F("hchi2", "track #chi^{2};track #chi^{2};# of events", 100, 0, 300);
  TH1F*  hchi2ndof = new TH1F("hchi2ndof", "track #chi^{2}/ndof;track #chi^{2}/n_{dof};# of events", 100, 0, 20);
  TH1F*  hchi2ndof_svt = new TH1F("hchi2_svt", "track #chi^{2}/dof svt;track #chi^{2}/dof;# of events", 100, 0, 20);
  TH1F*  hchi2ndof_bmtz = new TH1F("hchi2_bmtz", "track #chi^{2}/dof bmtz;track #chi^{2}/dof;# of events", 100, 0, 20);
  TH1F*  hchi2ndof_bmtc = new TH1F("hchi2_bmtc", "track #chi^{2}/dof bmtc;track #chi^{2}/dof;# of events", 100, 0, 20);
  TH1F*  hchi2prob = new TH1F("hchi2prob", "track #chi^{2} probability;p(#chi^2,n_{dof};# of events", 100, 0, 1);
  TH1F*  hndof = new TH1F("hndof", "degrees of freedom;n_{dof};# of events", 20, 0, 20);
  TH1F*  hncross = new TH1F("hncross", "number of crosses used;# of crosses;# of events", 20, 0, 20);
  TH1F*  hres = new TH1F("hres", "resolution;resolution (mm);# of clusters", 100, 0, 0.2);
  TH2F*  hresmodule =	new TH2F("hresmodule", "resolution (mm);resolution;module id", 100, 0, 0.2,103,0,103);
  TH1F*  hmodule = new TH1F("hmodule", "module;module id;# of crosses", 102, 0, 102);
  
  TH1F*  htrackderiv1 = new TH1F("htrackderiv1", "Track derivatives (doca);B_{i,0} (dimensionless);# of clusters", 100, -2, 2);
  TH1F*  htrackderiv2 = new TH1F("htrackderiv2", "Track derivatives (phi0);B_{i,1} (mm);# of clusters", 100, -300, 300);
  TH1F*  htrackderiv3 = new TH1F("htrackderiv3", "Track derivatives (z0);B_{i,2} (dimensionless);# of clusters", 100, -1.3, 1.3);
  TH1F*  htrackderiv4 = new TH1F("htrackderiv4", "Track derivatives (tan dep);B_{i,3} (mm);# of clusters", 100, -30, 30);
  TH2F*  htrackderiv1mod = new TH2F("htrackderiv1mod", "Track derivatives (doca);B_{i,0} (dimensionless);module id", 10, -4, 4,103,0,103);
  TH2F*  htrackderiv2mod = new TH2F("htrackderiv2mod", "Track derivatives (phi0);B_{i,1} (mm);module id", 100, -300, 300,103,0,103);
  TH2F*  htrackderiv3mod = new TH2F("htrackderiv3mod", "Track derivatives (z0);B_{i,2} (dimensionless);module id", 100, -1.3, 1.3,103,0,103);
  TH2F*  htrackderiv4mod = new TH2F("htrackderiv4mod", "Track derivatives (tan dep);B_{i,3} (mm);module id", 100, -30, 30,103,0,103);

  TH1F*  halignderiv1 = new TH1F("halignderiv1", "Alignment derivatives (shift x);A_{i,0} (dimensionless);# of clusters", 100, -4, 4);
  TH1F*  halignderiv2 = new TH1F("halignderiv2", "Alignment derivatives (shift y);A_{i,1} (dimensionless);# of clusters", 100, -4, 4);
  TH1F*  halignderiv3 = new TH1F("halignderiv3", "Alignment derivatives (shift z);A_{i,2} (dimensionless);# of clusters", 100, -0.03, 0.07);
  TH1F*  halignderiv4 = new TH1F("halignderiv4", "Alignment derivatives (rotate x);A_{i,3} (mm/rad);# of clusters", 100, -300, 300);
  TH1F*  halignderiv5 = new TH1F("halignderiv5", "Alignment derivatives (rotate y);A_{i,4} (mm/rad);# of clusters", 100, -300, 300);
  TH1F*  halignderiv6 = new TH1F("halignderiv6", "Alignment derivatives (rotate z);A_{i,5} (mm/rad);# of clusters", 100, -300, 300);

  TH2F*  halignderiv1mod = new TH2F("halignderiv1mod", "Alignment derivatives (shift x);A_{i,0} (dimensionless);module id", 100, -4, 4,103,0,103);
  TH2F*  halignderiv2mod = new TH2F("halignderiv2mod", "Alignment derivatives (shift y);A_{i,1} (dimensionless);module id", 100, -1, 1, 103,0,103);
  TH2F*  halignderiv3mod = new TH2F("halignderiv3mod", "Alignment derivatives (shift z);A_{i,2} (dimensionless);module id", 100, -0.03, 0.07,103,0,103);
  TH2F*  halignderiv4mod = new TH2F("halignderiv4mod", "Alignment derivatives (rotate x);A_{i,3} (mm/rad);module id", 100, -300, 300,103,0,103);
  TH2F*  halignderiv5mod = new TH2F("halignderiv5mod", "Alignment derivatives (rotate y);A_{i,4} (mm/rad);module id", 100, -300, 300,103,0,103);
  TH2F*  halignderiv6mod = new TH2F("halignderiv6mod", "Alignment derivatives (rotate z);A_{i,5} (mm/rad);module id", 100, -300, 300,103,0,103);

  TH1F*  hmeas = new TH1F("hmeas", "1D position (measured);1D position (mm);# of clusters", 100, -50, 50);
  TH1F*  hextrap = new TH1F("hextrap", "1D position (extrap);1D position (mm);# of clusters", 100, -50, 50);
  TH2F*  hmeasextrap = new TH2F("hmeasextrap", "1D position (extrap);meas. 1D position (mm);extrap. 1D position (mm)", 100, -50, 50,100,-50,50);

  TH1F*  hresidual = new TH1F("hresidual", "1D position residual;meas - extrap (mm);# of clusters", 100, -maxResid, maxResid);
  TH2F*  hresidualmod = new TH2F("hresidualmod", "1D position residual;meas - extrap (mm);module #", 100, -maxResid, maxResid,103,0,103);
  
  TH1F*  residuals_beamspot = new TH1F("hresidual_beamspot", "1D position residual (beamspot);meas - extrap (mm);# of clusters", 100, -maxResid, maxResid);
  
  vector<TH1*> residuals_BMT;
  for(int i = 0; i<18;i++){
    TH1F*  hresiduali = new TH1F(Form("hresidual_bmt_%d",i), "1D position residual BMT;meas - extrap (mm);# of clusters", 100, -maxResid, maxResid);
    residuals_BMT.push_back(hresiduali);
  }
  
  TH1F*  hresidual_SVT = new TH1F("hresidual_svt", "1D position residual (all SVT);meas - extrap (mm);# of clusters", 100, -maxResid, maxResid);
  
  
  
  TH1F*  hresidualnorm = new TH1F("hresidualnorm", "1D position residual normalized ;meas - extrap;# of clusters", 100, -3, 3);
  TH2F*  hresidualnormmod = new TH2F("hresidualnormmod", "1D position residual normalized;meas - extrap;module #", 100, -3, 3,103,0,103);

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
    double ndof_svt = 0;
    double ndof_bmtz = 0;
    double ndof_bmtc = 0;
    double chi2_svt = 0;
    double chi2_bmtz = 0;
    double chi2_bmtc = 0;
    
    for(int j = 0; j < aevent->GetMeasuredCovariance()->GetNrows(); j++){
      double res = sqrt((*aevent->GetMeasuredCovariance())[j][j]);
      double module = aevent->GetIndex()->At(j);
      hres->Fill(res);
      hmodule->Fill(module);
            
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
      double resid = (*aevent->GetMeasurements())(j)-(*aevent->GetTrackPrediction())(j);
      hresidual->Fill(resid);
      if(module >=84 && module < 102){
        if (module <=86 || module>=93 && module<=95 || module>=99){
          ndof_bmtc ++;
          chi2_bmtc+=resid*resid/(*aevent->GetMeasuredCovariance())[j][j];
        } else {
          ndof_bmtz ++;
          chi2_bmtz+=resid*resid/(*aevent->GetMeasuredCovariance())[j][j];
          
        }
        residuals_BMT[module-84]->Fill((*aevent->GetMeasurements())(j)-(*aevent->GetTrackPrediction())(j));
      } else if (module == 102){
        residuals_beamspot->Fill((*aevent->GetMeasurements())(j)-(*aevent->GetTrackPrediction())(j));
      } else {
        ndof_svt ++;
        chi2_svt += resid*resid/(*aevent->GetMeasuredCovariance())[j][j];
        hresidual_SVT->Fill((*aevent->GetMeasurements())(j)-(*aevent->GetTrackPrediction())(j));
      }
      hresidualmod->Fill((*aevent->GetMeasurements())(j)-(*aevent->GetTrackPrediction())(j),module);
      hresidualnorm->Fill(((*aevent->GetMeasurements())(j)-(*aevent->GetTrackPrediction())(j))/res);
      hresidualnormmod->Fill(((*aevent->GetMeasurements())(j)-(*aevent->GetTrackPrediction())(j))/res,module);
      
    }
    //hchi2ndof_svt->Fill(chi2_svt/ndof_svt);
    //hchi2ndof_bmtz->Fill(chi2_bmtz/ndof_bmtz);
    //hchi2ndof_bmtc->Fill(chi2_bmtc/ndof_bmtc);
    
    hchi2ndof_svt->Fill(chi2_svt/ndof);
    hchi2ndof_bmtz->Fill(chi2_bmtz/ndof);
    hchi2ndof_bmtc->Fill(chi2_bmtc/ndof);
  }
  
  //this is specific to the CVT
  double R = 150;
  int Nbins = 200;
  TH2F*  hmodulexy = createPlotXY(hmodule,"Modules", "clusters", "count");
    //new TH2F("hresmodulexy", "modules xy;x (mm);y (mm)", Nbins,-R,R,Nbins,-R,R);
  /*for(int m = 0; m<42; m++){
    
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
  */
  
  if(!(plotsDir==TString())){
    TLine* l= new TLine(); l->SetLineColor(kRed);
    TCanvas* c = new TCanvas("canvas","canvas",800,600);
    hchi2->Draw();drawLabel(label);c->SaveAs(plotsDir+"/chi2.png");
    hchi2ndof->Draw();drawLabel(label);c->SaveAs(plotsDir+"/chi2ndof.png");
    hchi2ndof_svt->Draw();drawLabel(label);c->SaveAs(plotsDir+"/chi2ndof_svt.png");
    hchi2ndof_bmtz->Draw();drawLabel(label);c->SaveAs(plotsDir+"/chi2ndof_bmtz.png");
    hchi2ndof_bmtc->Draw();drawLabel(label);c->SaveAs(plotsDir+"/chi2ndof_bmtc.png");
    hchi2prob->Draw();drawLabel(label);c->SaveAs(plotsDir+"/chi2prob.png");
    hndof->Draw();drawLabel(label);c->SaveAs(plotsDir+"/ndof.png");
    hncross->Draw();drawLabel(label);c->SaveAs(plotsDir+"/ncross.png");
    hres->Draw();drawLabel(label);c->SaveAs(plotsDir+"/resolution.png");
    hresmodule->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/res_vs_module.png");
    
    hmodule->Draw();

    l->DrawLine(10,0,10,hmodule->GetMaximum());
    l->DrawLine(24,0,24,hmodule->GetMaximum());
    l->DrawLine(42,0,42,hmodule->GetMaximum());
    TText * t = new TText();
    t->SetTextColor(kRed);
    t->DrawText(3,hmodule->GetMaximum()/10,"R1");
    t->DrawText(15,hmodule->GetMaximum()/10,"R2");
    t->DrawText(31,hmodule->GetMaximum()/10,"R3");
    drawLabel(label);c->SaveAs(plotsDir+"/module.png");
    
    htrackderiv1->Draw();drawLabel(label);c->SaveAs(plotsDir+"/B1.png");
    htrackderiv2->Draw();drawLabel(label);c->SaveAs(plotsDir+"/B2.png");
    htrackderiv3->Draw();drawLabel(label);c->SaveAs(plotsDir+"/B3.png");
    htrackderiv4->Draw();drawLabel(label);c->SaveAs(plotsDir+"/B4.png");
    htrackderiv1mod->SetStats(0);htrackderiv1mod->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/B1mod.png");
    htrackderiv2mod->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/B2mod.png");
    htrackderiv3mod->SetStats(0);htrackderiv3mod->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/B3mod.png");
    htrackderiv4mod->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/B4mod.png");
    gStyle->SetPadRightMargin(.2);
    createPlotXY(htrackderiv1mod,"B1 (doca)","B1 (dimensionless)")->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/B1xy.png");
    createPlotXY(htrackderiv2mod,"B2 (phi)","B2 (dimensionless)")->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/B2xy.png");
    createPlotXY(htrackderiv3mod,"B3 (z)","B3 (mm)")->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/B3xy.png");
    createPlotXY(htrackderiv4mod,"B4 (tandip)","B4 (mm)")->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/B4xy.png");
    gStyle->SetPadRightMargin(.10);
    halignderiv1->Draw();drawLabel(label);c->SaveAs(plotsDir+"/A1.png");
    halignderiv2->Draw();drawLabel(label);c->SaveAs(plotsDir+"/A2.png");
    halignderiv3->Draw();drawLabel(label);c->SaveAs(plotsDir+"/A3.png");
    halignderiv4->Draw();drawLabel(label);c->SaveAs(plotsDir+"/A4.png");
    halignderiv5->Draw();drawLabel(label);c->SaveAs(plotsDir+"/A5.png");
    halignderiv6->Draw();drawLabel(label);c->SaveAs(plotsDir+"/A6.png");
    halignderiv1mod->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/A1mod.png");
    halignderiv2mod->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/A2mod.png");
    halignderiv3mod->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/A3mod.png");
    halignderiv4mod->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/A4mod.png");
    halignderiv5mod->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/A5mod.png");
    halignderiv6mod->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/A6mod.png");

    gStyle->SetPadRightMargin(.2);
    createPlotXY(halignderiv1mod,"A1","A1 (dimensionless)")->Draw("COLZ1"); drawLabel(label);c->SaveAs(plotsDir+"/A1xy.png");
    createPlotXY(halignderiv2mod,"A2","A2 (dimensionless)")->Draw("COLZ1"); drawLabel(label);c->SaveAs(plotsDir+"/A2xy.png");
    createPlotXY(halignderiv3mod,"A3","A3 (dimensionless)")->Draw("COLZ1"); drawLabel(label);c->SaveAs(plotsDir+"/A3xy.png");
    createPlotXY(halignderiv4mod,"A4","A4 (mm)")->Draw("COLZ1"); drawLabel(label);c->SaveAs(plotsDir+"/A4xy.png");
    createPlotXY(halignderiv5mod,"A5","A5 (mm)")->Draw("COLZ1"); drawLabel(label);c->SaveAs(plotsDir+"/A5xy.png");
    createPlotXY(halignderiv6mod,"A6","A6 (mm)")->Draw("COLZ1"); drawLabel(label);c->SaveAs(plotsDir+"/A6xy.png");
    gStyle->SetPadRightMargin(.10);
    
    hmodulexy->Draw("COLZ1");drawLabel(label);c->SaveAs(plotsDir+"/modulexy.png");
    
    hmeas->SetTitle("1D hit positions");
    hmeas->Draw();hextrap->SetLineColor(kRed);hextrap->Draw("SAME");
    TLegend *legend = new TLegend();
    legend->AddEntry(hmeas, "measured");
    legend->AddEntry(hextrap, "track extrap");
    legend->Draw();
    drawLabel(label);c->SaveAs(plotsDir+"/meas_extrap.png");
    
    //BMT residuals plot
    for(int sec= 0;sec<3;sec++){
      int colors[] = {kRed,kCyan+1,kYellow+1,kGreen+1,kMagenta,kBlue};
      TLegend *legend = new TLegend(0.80, 0.6,1,0.95);
      int max = 0;
      for(int lay = 0; lay<6;lay++){
        auto h = residuals_BMT[lay*3+sec];
        h->SetLineColor(colors[lay]);
        h->SetTitle(Form("Residuals BMT Sec %d",sec));
        h->SetLineWidth(2);
        h->Draw(lay == 0 ? "" : "SAME");
        legend->AddEntry(h, Form("layer %d",lay));
        if(h->GetMaximum()>max)
          max = h->GetMaximum();
        /*TF1 *f1 = new TF1(Form("f%d",lay),"gaus",h->GetMean(),-5,5);
        f1->SetLineColor(colors[lay]);
        f1->SetLineStyle(7);
        f1->SetParameters(h->GetMaximum(), h->GetMean(), h->GetRMS() );*/
        h->Fit("gaus");
      }
      for(int lay = 0; lay<6;lay++){
        auto h = residuals_BMT[lay*3+sec];
        h->GetYaxis()->SetRangeUser(0,max);
      }
      legend->Draw();
      drawLabel(label);c->SaveAs(plotsDir+Form("/residuals_BMT_sec%d.png",sec));
    }
    
    //BMT residuals 18 panels
    


    residuals_beamspot->Draw();
    drawLabel(label);c->SaveAs(plotsDir+"/residuals_beamspot.png");
    hresidual_SVT->Draw();
    drawLabel(label);c->SaveAs(plotsDir+"/residuals_svt.png");
    
    hmeasextrap->Draw("COLZ1");l->DrawLine(-50,-50,50,50);drawLabel(label);c->SaveAs(plotsDir+"/meas_vs_extrap.png");
    hresidual->Draw();l->DrawLine(0,0,0,hresidual->GetMaximum());drawLabel(label);c->SaveAs(plotsDir+"/residual.png");
    
    hresidualmod->Draw("COLZ1");l->DrawLine(0,0,0,100);c->SetLogz(1);drawLabel(label);c->SaveAs(plotsDir+"/residualmod.png");
    c->SetLogz(0);
    hresidualnorm->Draw();l->DrawLine(0,0,0,hresidual->GetMaximum());drawLabel(label);c->SaveAs(plotsDir+"/residualnorm.png");
    hresidualnormmod->Draw("COLZ1");l->DrawLine(0,0,0,100);c->SetLogz(1);drawLabel(label);c->SaveAs(plotsDir+"/residualnormmod.png");
    c->SetLogz(0);
    
    TH1F*  hresidual_shift = new TH1F("hresidual_shift", "1D position residual;meas - extrap (mm);# of clusters", 100, -2, 2);
  
    
    gStyle->SetPadRightMargin(.2);
    createPlotXY(hresidualmod,"residuals","residual (mm)")->Draw("COLZ1"); drawLabel(label);c->SaveAs(plotsDir+"/residualxy.png");
    
  }
  
  
  //auto finish = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<double> elapsed = finish - start;
  //std::cout << "Elapsed time: " << elapsed.count()<< "s, events = "<<events<< "\n";
  
}
