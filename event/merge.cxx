#include <iostream>

#include "AlignEvent.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

int main(int argc, char * argv[]) 
{
  if(argc<2){
    cout << "usage:\nmerge out.root in1.root in2.root ..." << endl;
    return 0;
  }

  const int nEvents = 10;
   // initialize ROOT
  //  TROOT root("root", "root");
  std::cout << "hello world" <<std::endl;

  TFile * outfile = TFile::Open(argv[1], "RECREATE");

  TList* list = new TList();
  AlignInfo* alignInfo;
  for(int i = 2;i<argc;i++){
    TFile * infile = TFile::Open(argv[i]);
    list->Add(infile->Get("AlignTree")->Clone());
    alignInfo = (AlignInfo*)infile->Get("AlignInfo")->Clone();
  }

  outfile->cd();
  TTree * AlignTree = TTree::MergeTrees(list);
  AlignTree->Write();
  alignInfo->Write();
  outfile->Close();
  
}
