#include <iostream>

#include "AlignEvent.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

int main(int argc, char * argv[]) 
{
  if(argc<2){
    cout << "usage:\nevent fileName.root\nprints the matrices for the first 10 events in fileName.root" << endl;
    return 0;
  }

  const int nEvents = 10;
   // initialize ROOT
  TROOT root("root", "root");
  std::cout << "hello world" <<std::endl;
  
  // create event
  AlignEvent * event = new AlignEvent;

  
  // prepare reading
  TFile f1(argv[1], "READ");
  f1.Print();
  f1.ls();
  AlignInfo * ainfo = (AlignInfo *) f1.Get("AlignInfo");
  if (ainfo == 0) {
    return 2;
  }
  cout << "Alignment requested for " 
       << ainfo->GetNAlignables() << " alignables and "
       << ainfo->GetNParameters() << " parameters." << endl;
  TTree * t = (TTree *) f1.Get("AlignTree");
  if (t == 0)
    return 1;
  t->Print();
  TBranch *branch = t->GetBranch("AlignEvent");
  branch->SetAddress(&event);
  int nentries = t->GetEntries();
  cout << "tree has " << nentries << " entries." << endl;
  for (int j = 0; j < nEvents; j++) {
    cout << "now get entry " << endl;
    t->GetEntry(j);
    cout << "now print" << endl;
    cout << "measurements" << endl;
    event->GetMeasurements()->Print();
    cout << "covariance" << endl;
    event->GetMeasuredCovariance()->Print();
    cout << "track derivatives" << endl;
    event->GetTrackDerivatives()->Print();
    cout << "alignment derivatives" << endl;
    event->GetAlignmentDerivatives()->Print();
    cout << "track prediction" << endl;
    event->GetTrackPrediction()->Print();
  }
  f1.Close();
  
}
