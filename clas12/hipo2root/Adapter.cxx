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
#include <TStyle.h>
#include <TLine.h>
#include <TText.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <TLegend.h>
#include "../../event/AlignEvent.h"


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
  
  
  /////////////////////////////////////
  //ignore this just getting file name!
  TString inputFile;
  TString outputFile;
  
  
  vector<TString> inputFiles;
  for(Int_t i=1;i<argc;i++){
    TString opt=argv[i];
    if (opt.Contains("--out=")){
      outputFile = opt(6,opt.Sizeof());
    } else{
      inputFile = opt.ReplaceAll("--in=","");
      inputFiles.push_back(opt);
    }
  }
  //if there is no input file
  if(inputFile==TString())  {
    std::cout << " *** please provide a file name..." << std::endl;
    exit(0);
  }
  /////////////////////////////////////
  
  TFile * out = TFile::Open(outputFile, "RECREATE");
  
  
  AlignEvent *aevent = new AlignEvent();
  TTree *AlignTree = new TTree("AlignTree","Alignment data", 0);
  AlignTree->Branch("AlignEvent", & aevent, 64000, 0);
  
  //save the track parameters
  float q0=0;
  float q1=0;
  float q2=0;
  float q3=0;
  
  AlignTree->Branch("q0",&q0,"q0/F");
  AlignTree->Branch("q1",&q1,"q1/F");
  AlignTree->Branch("q2",&q2,"q2/F");
  AlignTree->Branch("q3",&q3,"q3/F");
  
  // auto save every MB
  AlignTree->SetAutoSave(1000000);
  
  
  // 42 SVT modules
  // TODO remove hardcode
  AlignInfo* ai = NULL;
  
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
    hipo::bank bank_q(factory.getSchema("Align::q"));
    hipo::bank bank_misc(factory.getSchema("Align::misc"));
    
    hipo::event event;
    while(r.next()){
      events++;
      r.read(event);
      event.getStructure(bank_misc);
      if(bank_misc.getRows() == 0)
        continue;
      event.getStructure(bank_B);
      event.getStructure(bank_V);
      event.getStructure(bank_m);
      event.getStructure(bank_c);
      event.getStructure(bank_I);
      event.getStructure(bank_A);
      event.getStructure(bank_q);
      for(int i = 0; i<bank_A.getRows(); i++){
        
        // 42 SVT modules
        // TODO remove hardcode
        if(ai==NULL){
          ai = new AlignInfo(bank_misc.getShort(bank_misc.getSchema().getEntryOrder("nalignables"),0),
                             bank_misc.getShort(bank_misc.getSchema().getEntryOrder("nparameters"),0));
          ai->Write();
          cout << "wrote align info" << endl;
        }
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
        //TVectorD* q;
        //getFromBank(q, bank_q,i);
        //cout << "checkpoint1" <<endl;
        cout << bank_q.getSchema().getEntryOrder("element_0_0") << " " << bank_q.getSchema().getEntryOrder("element_1_0") << " " << bank_q.getSchema().getEntryOrder("element_2_0") << " " << bank_q.getSchema().getEntryOrder("element_3_0") << endl;
        q0 = bank_q.getFloat(bank_q.getSchema().getEntryOrder("element_0_0"),i);
        q1 = bank_q.getFloat(bank_q.getSchema().getEntryOrder("element_1_0"),i);
        q2 = bank_q.getFloat(bank_q.getSchema().getEntryOrder("element_2_0"),i);
        q3 = bank_q.getFloat(bank_q.getSchema().getEntryOrder("element_3_0"),i);
        cout << q0 << " " << q1 << " " << q2 << " " << q3 << endl;
        //cout << "checkpoint2" <<endl;
        AlignTree->Fill();
        
        
      }
    }
    
    //r.close();
  }
  
  AlignTree->Write();
  
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Elapsed time: " << elapsed.count()<< "s, events = "<<events<< "\n";
  out->Write();
  out->Close();
}
