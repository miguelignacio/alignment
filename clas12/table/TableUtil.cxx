#include <cstdlib>
#include <list>
#include <iostream>
#include <stdio.h>
#include <iostream>
#include <fstream>
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
#include <string>


#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"

using namespace std;

int main(int argc, char * argv[]) {
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();
  
  TString inputFileName;
  
  TString oldFile;
  TString newFile;
  TString config;
  
  vector<TString> inputFiles;

  
  bool mergeTB = 0;
  for(Int_t i=1;i<argc;i++){
    TString opt=argv[i];
    if((opt.Contains("--in="))){
      inputFileName=opt(5,opt.Sizeof());  
    } 
    
    if (opt.Contains("--old=")){
      oldFile = opt(6,opt.Sizeof());
    }
    
    if (opt.Contains("--new=")){
      newFile = opt(6,opt.Sizeof());
    }

    if (opt.Contains("--config=")){
      config = opt(9,opt.Sizeof());
    }

    if (opt.Contains("--mergeTB")){
      mergeTB=1;
    }
  }
  //if there is no input file
  if(inputFileName==TString())  {
    std::cout << " *** please provide a file name..." << std::endl;
    exit(0);
  }
  /////////////////////////////////////

  //parse the config string
  
  int orderTx = -1;
  int orderTy = -1;
  int orderTz = -1;
  int orderRx = -1;
  int orderRy = -1;
  int orderRz = -1;

  int nparams = config.Sizeof()/2;
  for(int i = 0; i<config.Sizeof()/2; i++){
    if (config(2*i,2*i+2) == "Tx")
      orderTx = i;
    else if (config(2*i,2*i+2) == "Ty")
      orderTy =	i;
    else if (config(2*i,2*i+2) == "Tz")
      orderTz = i;
    else if (config(2*i,2*i+2) == "Rx")
      orderRx = i;
    else if (config(2*i,2*i+2) == "Ry")
      orderRy = i;
    else if (config(2*i,2*i+2) == "Rz")
      orderRz = i;
    else{
      cout << "invalid config string: "  << config << ".  Example:  TxTyRz -> translations in x, translations in y, rotation in z"<<  endl;
      return(0);
    }
  }
  
  
  
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

  //C.Print();
  TH1* hparams = new TH1F("params","params",N, 0, N);
  double params[N];
  double dparams[N];  
  for(int i = 0; i< N;i++){
    h=(TH1F*)inputFile->Get(Form("hAliPar%d",i));
    double n = h->GetBinContent(last);
    double dn = sqrt((*C)(i,i));
    cout << i <<" " << n<< " +- " << dn <<endl;
    hparams->SetBinContent(i+1,n);
    hparams->SetBinError(i+1,dn);
    params[i] = n;
    dparams[i] = dn;
  }
  
  double Tx[84];
  double Ty[84];
  double Tz[84];
  double Rx[84];
  double Ry[84];
  double Rz[84];
  double Ra[84];
  for(int i = 0; i<84; i++){
    if(orderTx >= 0)
      Tx[i] = params[nparams*i+orderTx];
    else
      Tx[i] = 0;
    if(orderTy >= 0)
      Ty[i] = params[nparams*i+orderTy];
    else
      Ty[i] = 0;
    if(orderTz >= 0)
      Tz[i] = params[nparams*i+orderTz];
    else
      Tz[i] = 0;
    if(orderRx >= 0)
      Rx[i] = params[nparams*i+orderRx];
    else
      Rx[i] = 0;
    if(orderRy >= 0)
      Ry[i] = params[nparams*i+orderRy];
    else
      Ry[i] = 0;
    if(orderRz >= 0)
      Rz[i] = params[nparams*i+orderRz];
    else
      Rz[i] = 0;
    Ra[i]=0;
  }
  
  ifstream fin (oldFile,ifstream::in);
  ofstream fout (newFile,ofstream::out);
  //char buffer
  TString s;
  
  double tx, ty, tz, rx, ry, rz, ra;

  if(mergeTB){
    for(int i = 0; i<42; i++){
      double * a[] = {Tx,Ty,Tz,Rx,Ry,Rz,Ra};
      for(int j = 0; j<7;j++){
	a[j][i]= (a[j][i]+a[j][i+42])/2;
	a[j][i+42] = a[j][i];
      }
    }
  }

  
  //three commented lines
  for(int i = 0; i<3; i++){
    string line;
    std::getline(fin, line);
    fout << line << endl;

  }
  
  int sector, layer, component;
  for(int i = 0; i<84; i++){
    fin >> sector >> layer >> component >> tx >> ty >> tz >> rx >> ry >> rz >> ra ;
    cout << "old tx = " << tx << ";  delta = " << Tx[i] << endl;
    int j = (layer % 2 ? 0 : 42) + sector-1;
    if(layer > 2)
      j += 10;
    if(layer > 4)
      j += 14;
    fout << sector << "\t" << layer << "\t" << component << "\t" << tx-Tx[j] << "\t" << ty-Ty[j] << "\t" << tz-Tz[j] << "\t" << rx-Rx[j] << "\t" << ry-Ry[j] << "\t" << rz-Rz[j] << "\t" << ra-Ra[j] << "\n";
    
  }
  fin.close();
  fout.close();
  cout << "wrote to " << newFile;
}
