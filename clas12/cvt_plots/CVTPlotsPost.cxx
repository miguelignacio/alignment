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
#include "TGaxis.h"

using namespace std;

void ReverseYAxis (TH1 *h)
{
   // Remove the current axis
   h->GetYaxis()->SetLabelOffset(999);
   h->GetYaxis()->SetTickLength(0);

   // Redraw the new axis
   gPad->Update();
   TGaxis *newaxis = new TGaxis(gPad->GetUxmin(),
                                gPad->GetUymax(),
                                gPad->GetUxmin()-0.001,
                                gPad->GetUymin(),
                                h->GetYaxis()->GetXmin(),
                                h->GetYaxis()->GetXmax(),
                                510,"+");
   newaxis->SetLabelOffset(-0.03);
  newaxis->SetLabelSize(h->GetYaxis()->GetLabelSize());
  newaxis->SetLabelFont(h->GetYaxis()->GetLabelFont());
   newaxis->Draw();
}

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
  gStyle->SetPalette(kViridis);
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();
  
  TString inputFileName;
  TString plotsDir;
  TString config="TxTy";
  
  vector<TString> inputFiles;
  for(Int_t i=1;i<argc;i++){
    TString opt=argv[i];
    if((opt.Contains("--in="))){
      inputFileName=opt(5,opt.Sizeof());

    } 
      
    if (opt.Contains("--plotsdir=")){
      plotsDir = opt(11,opt.Sizeof());
    }
    
    if (opt.Contains("--config=")){
      config = opt(9,opt.Sizeof());
    }
  }
  
  int orderTx = -1;
  int orderTy = -1;
  int orderTz = -1;
  int orderRx = -1;
  int orderRy = -1;
  int orderRz = -1;
  
  int nparams = config.Sizeof()/2;
  for(int i = 0; i<config.Sizeof()/2; i++){
    TString varname =config(2*i,2);
    if (varname == "Tx")
      orderTx = i;
    else if (varname == "Ty")
      orderTy =  i;
    else if (varname == "Tz")
      orderTz = i;
    else if (varname == "Rx")
      orderRx = i;
    else if (varname == "Ry")
      orderRy = i;
    else if (varname == "Rz")
      orderRz = i;
    else{
      cout << "invalid variable in config string: "<< varname << endl;
      cout << "invalid config string: "  << config << ".  Example:  TxTyRz -> translations in x, translations in y, rotation in z"<<  endl;
      return(0);
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
  TH1* hparamsz = new TH1F("paramsz","paramsz",102, 0, 102);
  TH1* hparamsrx = new TH1F("paramsrx","paramsrx",102, 0, 102);
  TH1* hparamsry = new TH1F("paramsry","paramsry",102, 0, 102);
  TH1* hparamsrz = new TH1F("paramsrz","paramsrz",102, 0, 102);
  
  //infer the number of parameters per module: if there is no BMT, then it's N/84, else it's N/102
  
  for(int i = 0; i< 102;i++){
    if(orderTx >= 0){
      hparamsx->SetBinContent(i+1,params[nparams*i+orderTx]);
      hparamsx->SetBinError(i+1,dparams[nparams*i+orderTx]);
    }
    if(orderTy >= 0){
      hparamsy->SetBinContent(i+1,params[nparams*i+orderTy]);
      hparamsy->SetBinError(i+1,dparams[nparams*i+orderTz]);
    }
    if(orderTz >= 0){
      hparamsz->SetBinContent(i+1,params[nparams*i+orderTz]);
      hparamsz->SetBinError(i+1,dparams[nparams*i]+orderTz);
    }
    if(orderRx >= 0){
      hparamsrx->SetBinContent(i+1,params[nparams*i+orderRx]);
      hparamsrx->SetBinError(i+1,dparams[nparams*i+orderRx]);
    }
    if(orderRy >= 0){
      hparamsry->SetBinContent(i+1,params[nparams*i+orderRy]);
      hparamsry->SetBinError(i+1,dparams[nparams*i+orderRy]);
    }
    if(orderRz >= 0){
      hparamsrz->SetBinContent(i+1,params[nparams*i+orderRz]);
      hparamsrz->SetBinError(i+1,dparams[nparams*i+orderRz]);
    }
  }
  if (orderTx >= 0){
    hparamsx->Draw();
    c->SaveAs(plotsDir+"/Tx.pdf");
    TH2F* hx_xy = createPlotXY(hparamsx, "translation x", "tx");
    hx_xy->Draw("COLZ");
    c->SaveAs(plotsDir+"/Tx_xy.pdf");
  }
  if (orderTy >= 0){
    hparamsy->Draw();
    c->SaveAs(plotsDir+"/Ty.pdf");
    TH2F* hy_xy = createPlotXY(hparamsy, "translation y", "ty");
    hy_xy->Draw("COLZ");
    c->SaveAs(plotsDir+"/Ty_xy.pdf");
  }
  if (orderTz >= 0){
    hparamsz->Draw();
    c->SaveAs(plotsDir+"/Tz.pdf");
    TH2F* hz_xy = createPlotXY(hparamsz, "translation z", "tz");
    hz_xy->Draw("COLZ");
    c->SaveAs(plotsDir+"/Tz_xy.pdf");
  }
  if (orderRx >= 0){
    hparamsrx->Draw();
    c->SaveAs(plotsDir+"/Rx.pdf");
    TH2F* hrx_xy = createPlotXY(hparamsrx, "rotation x", "rx");
    hrx_xy->Draw("COLZ");
    c->SaveAs(plotsDir+"/Rx_xy.pdf");
  }
  if (orderRy >= 0){
    hparamsry->Draw();
    c->SaveAs(plotsDir+"/Ry.pdf");
    TH2F* hry_xy = createPlotXY(hparamsry, "rotation y", "ry");
    hry_xy->Draw("COLZ");
    c->SaveAs(plotsDir+"/Ry_xy.pdf");
  }
  if (orderRz >= 0){
    hparamsrz->Draw();
    c->SaveAs(plotsDir+"/Rz.pdf");
    TH2F* hrz_xy = createPlotXY(hparamsrz, "rotation y", "rz");
    hrz_xy->Draw("COLZ");
    c->SaveAs(plotsDir+"/Rz_xy.pdf");
  }
  
  
  
  
  
  
  

  
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
  
  if(orderTx >=0 && orderTy>=0){
    double x[102];
    double y[102];
    double dx[102];
    double dy[102];
    for(int i = 0; i<102; i++){
      x[i] = params[nparams*i+orderTx];
      y[i] = params[nparams*i+orderTy];
      dx[i] = dparams[nparams*i+orderTx];
      dy[i] = dparams[nparams*i+orderTy];
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
  
  char* param_names[] = {"Tx","Ty","Tz","Rx","Ry","Rz"};
  for(int j = 0; j<nparams; j++){
    char* param_name = "";
    if(j == orderTx)
      param_name = "Tx";
    else if(j == orderTy)
      param_name = "Ty";
    else if(j == orderTz)
      param_name = "Tz";
    else if(j == orderRx)
      param_name = "Rx";
    else if(j == orderRy)
      param_name = "Ry";
    else if(j == orderRz)
      param_name = "Rz";
    
    for (int sec = 0; sec<3;sec++){
      TLegend *legend = new TLegend(0.75, 0.6,1,0.95);
      double max = 0;
      double min = 0;
      int colors[] = {kRed,kCyan+1,kYellow+1,kGreen+1,kMagenta,kBlue};
      for (int lay = 0; lay<6; lay++){
        int i = nparams*(3*lay+sec+84)+j;
        //cout << "nparam" << nparams << endl;
        TH1* h = (TH1*)inputFile->Get(Form("hAliPar%d",i));
        //cout << "h" << i << endl;
        
        h->SetTitle(Form("Evolution of %s parameters (BMT sec %d);track #;", param_name,sec));
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
        int i = nparams*(3*lay+sec+84)+j;
        auto h = (TH1*)inputFile->Get(Form("hAliPar%d",i));
        h->GetYaxis()->SetRangeUser(min,max);
        h->GetXaxis()->SetRangeUser(0,last);
      }
      legend->Draw();
      //save space by exporting these to png
      c->SaveAs(plotsDir+Form("/bmt_params_%s_sec%d.png", param_name,sec));
    }
  }
  
  
  int matsize = (18+84)*6;
  TH2* hcorr = new TH2D("correlations","correlation matrix;column;row", matsize,0,matsize,matsize,0,matsize);
  TH2* hcorr_bmt = new TH2D("correlations_bmt","correlations (bmt);column;row", matsize-84*nparams, 84*nparams,matsize,matsize-84*nparams,84*nparams,matsize);
  hcorr->SetStats(0);
  hcorr->GetZaxis()->SetLabelSize(0.022);
  hcorr_bmt->SetStats(0);
  hcorr_bmt->GetZaxis()->SetLabelSize(0.022);
  hcorr->SetMaximum(1);
  hcorr->SetMinimum(-1);
  hcorr_bmt->SetMaximum(1);
  hcorr_bmt->SetMinimum(-1);
  auto cc = *C;
  for(int i = 0; i<matsize; i++){
    double cii =cc(i,i);
    for(int j = 0; j<matsize; j++){
      double cjj = cc(j,j);
      auto corr = cc(i,j)/sqrt(cii*cjj);
      hcorr->SetBinContent(1+i,matsize-j,corr);
      if(i>=84*nparams)
        hcorr_bmt->SetBinContent(i-84*nparams+1,matsize-j,corr);
    }
                         
  }
  
  
  TText* t = new TText();
  t->SetTextAlign(22);
  t->SetTextColor(kRed);
  t->SetTextColorAlpha(kRed,0.5);
  
  // create a color scheme in which black is zero, magenta is positive,
  // cyan is negative
  /*int n = 3;
  double r[3] = {0,0,1};
  double g[3] = {1,0,0};
  double b[3] = {1,0,1};
  double s[3] = {0,0.5,1.0};*/
  int n = 5;
  double r[5] = {0.4,0,0,1,1};
  double g[5] = {1,1,0,0,0.4};
  double b[5] = {1,1,0,1,1};
  double s[5] = {0,0.13,0.5,0.87,1.0};
  TColor::CreateGradientColorTable(n, s, r, g, b, 99);
  //c->SetLogz();
  hcorr->Draw("COLZ");
  t->SetTextSize(.15);
  t->DrawText(N*42./102,matsize-N*42./102,"SVT");
  t->SetTextSize(.05);
  t->DrawText(N*93./102,matsize-N*93./102,"BMT");
  TLine * l = new TLine();
  l->SetLineColorAlpha(kRed,1);
  l->SetLineStyle(3);
  l->DrawLine(0, matsize-matsize*84./102, matsize, matsize-matsize*84./102);
  l->DrawLine(matsize*84./102, 0, matsize*84./102, matsize);
  
  //hopefully this makes each bin exactly 2 pixels
  c->SetRightMargin(0.10);
  c->SetLeftMargin(0.10);
  c->SetTopMargin(0.10);
  c->SetBottomMargin(0.10);
  c->SetCanvasSize(matsize*10/4,matsize*10/4);
  ReverseYAxis(hcorr);
  c->SaveAs(plotsDir+"/correlations.png");
  
  hcorr_bmt->Draw("COLZ");
  ReverseYAxis(hcorr_bmt);
  //hopefully this makes each bin exactly 2 pixels
  c->SetRightMargin(0.10);
  c->SetLeftMargin(0.10);
  c->SetTopMargin(0.10);
  c->SetBottomMargin(0.10);
  c->SetCanvasSize(N*10/4,N*10/4);
  c->SaveAs(plotsDir+"/correlations_bmt.png");
  
  /*TLegend* legend = new TLegend(0.75, 0.6,1,0.95);
  
  int colors[] = {kRed,kCyan+1,kYellow+1,kGreen+1,kMagenta,kBlue};
  for (int j = 0; j<6; j++){
    int i =j+168+36;
    TH1* h = (TH1*)inputFile->Get(Form("hAliPar%d",i));
    h->SetTitle("Evolution of parameters (beamspot);track #");
    h->SetLineColor(colors[j]);
    h->SetLineWidth(2);
    h->Draw(j == 0 ?"":"SAME");
    legend->AddEntry(h, param_names[j]);
  }
  
  for (int j = 0; j<6; j++){
    int i =j+168+36;
    TH1* h = (TH1*)inputFile->Get(Form("hAliPar%d",i));
    h->GetXaxis()->SetRangeUser(0,last);
    h->GetYaxis()->SetRangeUser(-4,3);
  }
    
  legend->Draw();
  c->SaveAs(plotsDir+"/params_beam_xy.pdf");
  */
  
  
  if(nparams==2){
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
