#include <cstdlib>
#include <list>
#include <iostream>
#include <chrono>
#include "hipo4/reader.h"
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <TLegend.h>
#include "../event/AlignEvent.h"


#include <iostream>



#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

void getFromBank(TMatrixD* mat, hipo::bank bank, int i){
  int rows = bank.getInt(bank.getSchema().getEntryOrder("rows"),i);
  int columns = bank.getInt(bank.getSchema().getEntryOrder("columns"),i);
  mat->ResizeTo(rows,columns);
  int eo = 0;
  //cout << bank.getAddress() << endl;
  for(int j = 0; j < rows; j++){
    for(int k = 0; k <columns; k++){
      double val = bank.getFloat(bank.getSchema().getEntryOrder(Form("element_%d_%d",j,k)),i);
      //cout <<" " << Form("%.7f",val);
      (*mat)(j,k) = val;
      eo++;
    }
    //cout << endl<<endl;
  }
}

void getFromBank(TMatrixDSym* mat, hipo::bank bank, int i){
  int rows = bank.getInt(bank.getSchema().getEntryOrder("rows"),i);
  int columns = bank.getInt(bank.getSchema().getEntryOrder("columns"),i);
  mat->ResizeTo(rows,columns);
  int eo = 0;
  //cout << bank.getAddress() << endl;
  for(int j = 0; j < rows; j++){
    for(int k = 0; k <columns; k++){
      double val = bank.getFloat(bank.getSchema().getEntryOrder(Form("element_%d_%d",j,k)),i);
      //cout <<" " << Form("%.7f",val);
      (*mat)(j,k) = val;
      eo++;
    }
    //cout << endl<<endl;
  }
}



void getFromBank(TArrayI* mat, hipo::bank bank, int i){
  int rows = bank.getInt(bank.getSchema().getEntryOrder("rows"),i);
  int columns = bank.getInt(bank.getSchema().getEntryOrder("columns"),i);
  mat->Set(rows);
  int eo = 0;
  //cout << "rows = "<< rows << endl;
  for(int j = 0; j < rows; j++){
    for(int k = 0; k <1; k++){
      double val = bank.getFloat(bank.getSchema().getEntryOrder(Form("element_%d_%d",j,k)),i);
      //cout <<" " << Form("%.7f",val);
      mat->SetAt(val,j);
      eo++;
    }
    //cout << endl<<endl;
  }
}

void getFromBank(TVectorD* mat, hipo::bank bank, int i){
  int rows = bank.getInt(bank.getSchema().getEntryOrder("rows"),i);
  int columns = bank.getInt(bank.getSchema().getEntryOrder("columns"),i);
  mat->ResizeTo(rows);
  int eo = 0;
  //cout << bank.getAddress() << endl;
  for(int j = 0; j < rows; j++){
    for(int k = 0; k <1; k++){
      double val = bank.getFloat(bank.getSchema().getEntryOrder(Form("element_%d_%d",j,k)),i);
      //cout <<" " << Form("%.7f",val);
      (*mat)(j) = val;
      eo++;
    }
    //cout << endl<<endl;
  }
}



int main(int argc, char * argv[]) {
  // Record start time
  auto start = std::chrono::high_resolution_clock::now();
  
  bool debug = 0;
  /////////////////////////////////////
  //ignore this just getting file name!
  TString inputFile;
  TString outputFile;
  
  vector<TString> inputFiles;
  for(Int_t i=1;i<argc;i++){
    TString opt=argv[i];
    if((opt.Contains("--in="))){
      inputFile=opt(5,opt.Sizeof());
      inputFiles.push_back(inputFile);
    } else if (opt.Contains("--out=")){
      outputFile = opt(6,opt.Sizeof());
    } if (opt.Contains("--debug")){
      debug = 1;
    }
  }
  //if there is no input file
  if(inputFile==TString())  {
    std::cout << " *** please provide a file name..." << std::endl;
    exit(0);
  }
  /////////////////////////////////////
  
  TFile * out = TFile::Open(outputFile, "RECREATE");
  
  // 42 SVT modules
  // TODO remove hardcode
  AlignInfo* ai = new AlignInfo(42,6);

  ai->Write();
  AlignEvent *aevent = new AlignEvent();
  TTree *AlignTree = new TTree("AlignTree","Alignment data", 0);
  AlignTree->Branch("AlignEvent", & aevent, 64000, 0);
  // auto save every MB
  AlignTree->SetAutoSave(1000000);

  // these histograms are for debugging.  The limits and labels for
  // some of these histograms are taylored for the Clas12 CVT.
  
  TH1F*  hchi2 = new TH1F("hchi2", "track #chi^{2}", 100, 0, 100);
  TH1F*  hchi2ndof = new TH1F("hchi2ndof", "track #chi^{2}/ndof;track #chi^{2}/n_{dof};# of events", 100, 0, 10);
  TH1F*  hchi2prob = new TH1F("hchi2prob", "track #chi^{2} probability;p(#chi^2,n_{dof};# of events", 100, 0, 1);
  TH1F*  hndof = new TH1F("hndof", "degrees of freedom;n_{dof};# of events", 20, 0, 20);
  TH1F*  hncross = new TH1F("hncross", "number of crosses used;# of crosses;# of events", 20, 0, 20);
  TH1F*  hres = new TH1F("hres", "resolution;resolution (mm);# of clusters", 100, 0, 0.1);
  TH2F*  hresmodule =	new TH2F("hresmodule", "resolution (mm);resolution;module id", 100, 0, 0.1,50,0,50);
  TH1F*  hmodule = new TH1F("hmodule", "module;module number i;# of clusters", 50, 0, 50);

  TH1F*  htrackderiv1 = new TH1F("htrackderiv1", "Track derivatives (x);B_{i,0};# of clusters", 100, -4, 4);
  TH1F*  htrackderiv2 = new TH1F("htrackderiv2", "Track derivatives (z);B_{i,1};# of clusters", 100, -0.07, 0.03);
  TH1F*  htrackderiv3 = new TH1F("htrackderiv3", "Track derivatives (tx);B_{i,2};# of clusters", 100, -300, 300);
  TH1F*  htrackderiv4 = new TH1F("htrackderiv4", "Track derivatives (tz);B_{i,3};# of clusters", 100, -10, 10);
  TH2F*  htrackderiv1mod = new TH2F("htrackderiv1mod", "Track derivatives (x);B_{i,0};module id", 100, -4, 4,50,0,50);
  TH2F*  htrackderiv2mod = new TH2F("htrackderiv2mod", "Track derivatives (z);B_{i,1};module id", 100, -0.07, 0.03,50,0,50);
  TH2F*  htrackderiv3mod = new TH2F("htrackderiv3mod", "Track derivatives (tx);B_{i,2};module id", 100, -300, 300,50,0,50);
  TH2F*  htrackderiv4mod = new TH2F("htrackderiv4mod", "Track derivatives (tz);B_{i,3};module id", 100, -10, 10,50,0,50);
  
  int events = 0, tracks=0;
   for(Int_t filenum=0;filenum<inputFiles.size();filenum++){
     //create the event reader
     inputFile = inputFiles[filenum];
     cout<<"Analysing hipo file "<<inputFile<<endl;
     hipo::reader r;
     r.open(inputFile);
     hipo::dictionary factory;
     r.readDictionary(factory);
     hipo::bank bank_A(factory.getSchema("Align::A"));
     hipo::bank bank_B(factory.getSchema("Align::B"));
     hipo::bank bank_V(factory.getSchema("Align::V"));
     hipo::bank bank_m(factory.getSchema("Align::m"));
     hipo::bank bank_c(factory.getSchema("Align::c"));
     hipo::bank bank_I(factory.getSchema("Align::I"));
     hipo::bank bank_misc(factory.getSchema("Align::misc"));
     
     hipo::event event;

     while(r.next()){
       events++;
       r.read(event);
       event.getStructure(bank_A);
       event.getStructure(bank_B);
       event.getStructure(bank_V);
       event.getStructure(bank_m);
       event.getStructure(bank_c);
       event.getStructure(bank_I);
       event.getStructure(bank_misc);
       
       for(int i = 0; i<bank_A.getRows(); i++){

	 //cout << bank_I.getFloat(2,0) << endl;
	 getFromBank(aevent->GetAlignmentDerivatives(),bank_A, i);
	 getFromBank(aevent->GetTrackDerivatives(),bank_B, i);
	 getFromBank(aevent->GetMeasurements(),bank_m, i);
         getFromBank(aevent->GetTrackPrediction(),bank_c, i);
         getFromBank(aevent->GetMeasuredCovariance(),bank_V, i);
         getFromBank(aevent->GetIndex(),bank_I, i);

	 int ndof = bank_misc.getShort(bank_misc.getSchema().getEntryOrder("ndof"),i);
	 float chi2 = bank_misc.getFloat(bank_misc.getSchema().getEntryOrder("chi2"),i); 
	 aevent->SetRun(bank_misc.getInt(bank_misc.getSchema().getEntryOrder("run"),i));
	 aevent->SetEvent(bank_misc.getInt(bank_misc.getSchema().getEntryOrder("event"),i));
	 aevent->SetChi2(bank_misc.getFloat(bank_misc.getSchema().getEntryOrder("chi2"),i));;
	 aevent->SetNdof(bank_misc.getShort(bank_misc.getSchema().getEntryOrder("ndof"),i));
	 aevent->SetTrackNumber(bank_misc.getInt(bank_misc.getSchema().getEntryOrder("track"),i));
	 AlignTree->Fill();

	 hchi2->Fill(chi2);
	 hchi2ndof->Fill(chi2/ndof);
	 hchi2prob->Fill(TMath::Prob(chi2, ndof));
	 hndof->Fill(ndof);
	 hncross->Fill(aevent->GetMeasuredCovariance()->GetNrows()/2);
	 for(int j = 0; j < aevent->GetMeasuredCovariance()->GetNrows(); j++){
	   double res = sqrt((*aevent->GetMeasuredCovariance())[j][j]);
	   int module = aevent->GetIndex()->At(j);
	   hres->Fill(res);
	   hresmodule->Fill(res,module);
	   hmodule->Fill(module);
	   TMatrixD *B = aevent->GetTrackDerivatives();
	   htrackderiv1->Fill((*B)[j][0]);
           htrackderiv2->Fill((*B)[j][1]);
           htrackderiv3->Fill((*B)[j][2]);
           htrackderiv4->Fill((*B)[j][3]);
	   htrackderiv1mod->Fill((*B)[j][0],module);
           htrackderiv2mod->Fill((*B)[j][1],module);
           htrackderiv3mod->Fill((*B)[j][2],module);
           htrackderiv4mod->Fill((*B)[j][3],module);
	 }
       }
     }
     
    
   }



   
   AlignTree->Write();

   hchi2->Write();
   hchi2ndof->Write();
   hchi2prob->Write();
   hndof->Write();
   hncross->Write();
   hres->Write();
   hresmodule->Write();
   hmodule->Write();

   htrackderiv1->Write();
   htrackderiv2->Write();
   htrackderiv3->Write();
   htrackderiv4->Write();
   htrackderiv1mod->Write();
   htrackderiv2mod->Write();
   htrackderiv3mod->Write();
   htrackderiv4mod->Write();
   
   auto finish = std::chrono::high_resolution_clock::now();
   std::chrono::duration<double> elapsed = finish - start;
   std::cout << "Elapsed time: " << elapsed.count()<< "s, events = "<<events<< "\n";
   out->Write();
   out->Close();
}
