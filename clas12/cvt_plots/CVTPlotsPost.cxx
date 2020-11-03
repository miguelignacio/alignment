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
#include <TMatrixTSym.h>
#include <TGraphErrors.h>
//#include "../event/AlignEvent.h"


#include <iostream>



#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"

using namespace std;


TH2F* createPlotXY(TH1* h, TString title, TString ztitle){
  
  double R = 150;
  int Nbins = 200;
  TH2F*  hxy =   new TH2F(h->GetName()+(TString)"_xy", title +";x (mm);y (mm);"+ztitle, Nbins,-R,R,Nbins,-R,R);

  
  //cout << prof->GetNbinsX() << "bins" << endl;
  for(int mm = 0; mm<92; mm++){
    int m = mm%50;
    double n = h->GetBinContent(mm+1);
    
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
  
  int last = 0;
  TH1* h = (TH1*)inputFile->Get(Form("hAliPar2"));
  int end = (int)h->GetXaxis()->GetXmax();
  for(int j=end; j>=0;j--) {
    if(h->GetBinContent(j)!=0){
      last=j;
      break;
    }
  }
  cout << "last event is " << last << endl;
  TH1* hparams = new TH1F("params","params",200, 0, 200);
  double params[168];
  double dparams[168];
  //auto C = *(TMatrixTSym*)inputFile->Get("CovarianceMatrix");
  TMatrixTSym<double> * C = (TMatrixTSym<double>*)inputFile->Get("CovarianceMatrix");
  //C.Print();
  for(int i = 0; i< 168;i++){
    h=(TH1F*)inputFile->Get(Form("hAliPar%d",i));
    double n = h->GetBinContent(last);
    double dn = sqrt((*C)(i,i));
    cout << i <<" " << n<< " +- " << dn <<endl;
    hparams->SetBinContent(i+1,n);
    hparams->SetBinError(i+1,dn);
    params[i] = n;
    dparams[i] = dn;
  }
  TCanvas* c = new TCanvas("canvas","canvas",800,600);
  hparams->Draw();
  c->SaveAs("params.pdf");
  TH1* hparamsx = new TH1F("paramsx","paramsx",100, 0, 100);
  TH1* hparamsy = new TH1F("paramsy","paramsy",100, 0, 100);
  for(int i = 0; i< 84;i++){
    hparamsx->SetBinContent(i+1+8*(i>42),params[2*i]);
    hparamsx->SetBinError(i+1+8*(i>42),dparams[2*i]);

    hparamsy->SetBinContent(i+1+8*(i>42),params[2*i+1]);
    hparamsy->SetBinError(i+1+8*(i>42),dparams[2*i+1]);
  }
  hparamsx->Draw();
  c->SaveAs(plotsDir+"/Tx.pdf");
  hparamsy->Draw();
  c->SaveAs(plotsDir+"/Ty.pdf");

  TH2F* hx_xy = createPlotXY(hparamsx, "translation x", "x");
  hx_xy->Draw("COLZ");
  c->SaveAs(plotsDir+"/Tx_xy.pdf");
  TH2F* hy_xy = createPlotXY(hparamsy, "translation y", "y");
  hy_xy->Draw("COLZ");
  c->SaveAs(plotsDir+"/Ty_xy.pdf");

  TH1F* allVals = new TH1F("all params", "all params",50, -10,10);
  for(int i = 0; i<168;i++){
    allVals->Fill(params[i]);
  }
  allVals->Draw();
  c->SaveAs(plotsDir+"/allVals.pdf");

  TH1F* diagCovs = new TH1F("diagonal cov elements", "sqrt(Cii)",50,-0.1,0.1);
  for(int i = 0; i<168;i++){
    diagCovs->Fill(dparams[i]);
  }
  diagCovs->Draw();
  c->SaveAs(plotsDir+"/diagonals.pdf");

  TH1F* offDiagCovs = new TH1F("off-diagonal cov elements", "Cij",50,-0.003,0.003);
  for(int i = 0; i<168;i++){
    for(int j = 0; j<i; j++){
      offDiagCovs->Fill((*C)(i,j));
    }
  }
  offDiagCovs->Draw();
  c->SaveAs(plotsDir+"/offDiagonals.pdf");

  TH1F* offDiagCovsNorm = new TH1F("off-diagonal cov elements normalized", "Cij/sqrt(Cii*Cjj)",50,-1,1);
  for(int i = 0; i<168;i++){
    for(int j = 0; j<i; j++){
      offDiagCovsNorm->Fill((*C)(i,j)/sqrt((*C)(i,i)*(*C)(j,j)));
    }
  }
  offDiagCovsNorm->Draw();
  c->SaveAs(plotsDir+"/offDiagonalsNormed.pdf");
  
  
  double x[84];
  double y[84];
  double dx[84];
  double dy[84];
  for(int i = 0; i<84; i++){
    x[i] = params[2*i];
    y[i] = params[2*i+1];
    dx[i] = dparams[2*i];
    dy[i] = dparams[2*i+1];
  }

  TGraphErrors * g =new TGraphErrors(84,x,y,dx,dy);
  g->Draw("AP");
  
  c->SaveAs(plotsDir+"/xvsy.pdf");



}
