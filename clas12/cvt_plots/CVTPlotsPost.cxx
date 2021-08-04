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

double angle(double phi){
  double pi = TMath::Pi();
  while(phi>pi)
    phi-=2*pi;
  while(phi<-pi)
    phi+=2*pi;
  return phi;
}

TH2F* createPlotXY(TH1* h, TString title, TString ztitle){
  
  double R = 235;
  int Nbins = 235;
  TH2F*  hxy =   new TH2F(h->GetName()+(TString)"_xy", title +";-x (mm);y (mm);"+ztitle, Nbins,-R,R,Nbins,-R,R);

  
  //cout << prof->GetNbinsX() << "bins" << endl;
  for(int mm = 0; mm<102; mm++){
    double n = h->GetBinContent(mm+1);
    if(mm < 84){
      int m = mm%42;
      
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
      cout<<  mm << " " << layer << " " << sector << " " << n << " " << phisector*180/TMath::Pi() << endl;

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

  //auto C = *(TMatrixTSym*)inputFile->Get("CovarianceMatrix");
  TMatrixTSym<double> * C = (TMatrixTSym<double>*)inputFile->Get("CovarianceMatrix");
  int N = C->GetNrows();
  cout << N << " parameters total" << endl;
  //C.Print();
  TH1* hparams = new TH1F("params","params",N, 0, N);
  double params[102*6];
  double dparams[102*6];  
  for(int i = 0; i< 102*6;i++){
    if(i<N){
      h=(TH1F*)inputFile->Get(Form("hAliPar%d",i));
      double n = h->GetBinContent(last);
      double dn = sqrt((*C)(i,i));
      cout << i <<" " << n<< " +- " << dn <<endl;
      hparams->SetBinContent(i+1,n);
      hparams->SetBinError(i+1,dn);
      params[i] = n;
      dparams[i] = dn;
    } else{
      params[i] = 0;
      dparams[i] = 0;
    }
  }
  TCanvas* c = new TCanvas("canvas","canvas",800,600);
  c->SetRightMargin(0.15);
  hparams->Draw();
  c->SaveAs("params.pdf");
  TH1* hparamsx = new TH1F("paramsx","paramsx",102, 0, 102);
  TH1* hparamsy = new TH1F("paramsy","paramsy",102, 0, 102);
  TH1* hparamsr = new TH1F("paramsr","paramsr",102, 0, 102);
  
  //infer the number of parameters per module: if there is no BMT, then it's N/84, else it's N/102
  int nparam = N%102 == 0? N/102 : N/84;
  for(int i = 0; i< 102;i++){
    hparamsx->SetBinContent(i+1,params[nparam*i]);
    hparamsx->SetBinError(i+1,dparams[nparam*i]);

    hparamsy->SetBinContent(i+1,params[nparam*i+1]);
    hparamsy->SetBinError(i+1,dparams[nparam*i+1]);
    if(nparam>2){
      hparamsr->SetBinContent(i+1,params[nparam*i+2]);
      hparamsr->SetBinError(i+1,dparams[nparam*i+2]);
    }
  }
  hparamsx->Draw();
  c->SaveAs(plotsDir+"/Tx.pdf");
  hparamsy->Draw();
  c->SaveAs(plotsDir+"/Ty.pdf");
  if(nparam>2){
    hparamsr->Draw();
    c->SaveAs(plotsDir+"/Rz.pdf");
  }
  TH2F* hx_xy = createPlotXY(hparamsx, "translation x", "x");
  hx_xy->Draw("COLZ");
  c->SaveAs(plotsDir+"/Tx_xy.pdf");
  TH2F* hy_xy = createPlotXY(hparamsy, "translation y", "y");
  hy_xy->Draw("COLZ");
  c->SaveAs(plotsDir+"/Ty_xy.pdf");

  TH2F* hr_xy = createPlotXY(hparamsr, "rotation z", "phi");
  hr_xy->Draw("COLZ");
  c->SaveAs(plotsDir+"/Rz_xy.pdf");

  
  TH1F* allVals = new TH1F("all params", "all params",50, -10,10);
  for(int i = 0; i<N;i++){
    allVals->Fill(params[i]);
  }
  allVals->Draw();
  c->SaveAs(plotsDir+"/allVals.pdf");

  TH1F* diagCovs = new TH1F("diagonal cov elements", "sqrt(Cii)",50,-0.1,0.1);
  for(int i = 0; i<N;i++){
    diagCovs->Fill(dparams[i]);
  }
  diagCovs->Draw();
  c->SaveAs(plotsDir+"/diagonals.pdf");

  TH1F* offDiagCovs = new TH1F("off-diagonal cov elements", "Cij",50,-0.003,0.003);
  for(int i = 0; i<N;i++){
    for(int j = 0; j<i; j++){
      offDiagCovs->Fill((*C)(i,j));
    }
  }
  offDiagCovs->Draw();
  c->SaveAs(plotsDir+"/offDiagonals.pdf");

  TH1F* offDiagCovsNorm = new TH1F("off-diagonal cov elements normalized", "Cij/sqrt(Cii*Cjj)",50,-1,1);
  for(int i = 0; i<N;i++){
    for(int j = 0; j<i; j++){
      offDiagCovsNorm->Fill((*C)(i,j)/sqrt((*C)(i,i)*(*C)(j,j)));
    }
  }
  offDiagCovsNorm->Draw();
  c->SaveAs(plotsDir+"/offDiagonalsNormed.pdf");
  
  {
    double x[102];
    double y[102];
    double dx[102];
    double dy[102];
    for(int i = 0; i<102; i++){
      x[i] = params[nparam*i];
      y[i] = params[nparam*i+1];
      dx[i] = dparams[nparam*i];
      dy[i] = dparams[nparam*i+1];
    }
    
    TGraphErrors * g =new TGraphErrors(102,x,y,dx,dy);
    g->Draw("AP"); 
    c->SaveAs(plotsDir+"/xvsy.pdf");
  }

  TH1* hstability = new TH1F("stability","stability;param # ;#Delta{p} in last 1/3 of data",N, 0, N);
  for(int i = 0; i< N;i++){
    h=(TH1F*)inputFile->Get(Form("hAliPar%d",i));
    double nlast = h->GetBinContent(last);
    double n23 = h->GetBinContent(last*2/3);
    hstability->SetBinContent(i+1,nlast-n23);
  }

  hstability->Draw();
  c->SaveAs(plotsDir+"/stability.pdf");
  

  for(int j = 0; j<2; j++){
    for (int sec = 0; sec<3;sec++){
      TLegend *legend = new TLegend(0.75, 0.6,1,0.95);
      double max = 0;
      double min = 0;
      int colors[] = {kRed,kCyan+1,kYellow+1,kGreen+1,kMagenta,kBlue};
      for (int lay = 0; lay<6; lay++){
        int i = 6*lay+2*sec+j+168;
        TH1* h = (TH1*)inputFile->Get(Form("hAliPar%d",i));
        h->SetTitle(Form("Evolution of %s parameters (BMT sec %d);track #;", j==1?"y":"x",sec));
        h->Draw(lay == 0 ? "": "SAME");
        h->SetLineWidth(2);
        h->SetLineColor(colors[lay]);
        legend->AddEntry(h, Form("layer %d",lay));
        if(h->GetMaximum()>max)
          max = h->GetMaximum();
        if(h->GetMinimum()<min)
          min = h->GetMinimum();
      }
      for(int lay = 0; lay<6;lay++){
        int i = 6*lay+2*sec+j+168;
        auto h = (TH1*)inputFile->Get(Form("hAliPar%d",i));
        h->GetYaxis()->SetRangeUser(min,max);
        h->GetXaxis()->SetRangeUser(0,last);
      }
      legend->Draw();
      c->SaveAs(plotsDir+Form("/bmt_params_%s_sec%d.pdf", j==1?"y":"x",sec));
    }
  }
  
  TLegend* legend = new TLegend(0.75, 0.6,1,0.95);
  
  int colors[] = {kRed,kCyan+1,kYellow+1,kGreen+1,kMagenta,kBlue};
  for (int j = 0; j<2; j++){
    int i =j+168+36;
    TH1* h = (TH1*)inputFile->Get(Form("hAliPar%d",i));
    h->SetTitle("Evolution of parameters (beamspot);track #");
    h->SetLineColor(colors[j]);
    h->SetLineWidth(2);
    h->Draw(j == 0 ?"":"SAME");
    legend->AddEntry(h, j == 0? "x":"y");
    
  }
  for (int j = 0; j<2; j++){
    int i =j+168+36;
    TH1* h = (TH1*)inputFile->Get(Form("hAliPar%d",i));
    h->GetXaxis()->SetRangeUser(0,last);
    h->GetYaxis()->SetRangeUser(-4,3);
  }
    
  legend->Draw();
  c->SaveAs(plotsDir+"/params_beam_xy.pdf");
  
  
  
  if(nparam==2){
    cout << "double tx = {" ;
    for(int i =0;i<N/2;i++){
      cout << params[2*i];
      if(i != N/2-1){
	cout << ", ";
	if(i%5 == 0)
	  cout << endl;
      }
      else cout << "};\n";
    }
    cout << "double ty = {" ;
    for(int i =0;i<N/2;i++){
      cout << params[2*i+1];
      if(i != N/2-1){
        cout << ", ";
        if(i%5 == 0)
          cout << endl;
      }
      else cout << "};\n";
    }
    
  }

}
