#include <iostream>

#include "AlignEvent.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

int main(int argc, char * argv[]) 
{
  const int nEvents = 10;
  /* // initialize ROOT
  TROOT root("root", "root");

  // create event
  AlignEvent * event = new AlignEvent;

  // prepare reading
  TFile f1("alignment_data.root", "READ");
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
    event->GetMeasurements()->Print();
    event->GetMeasuredCovariance()->Print();
    event->GetTrackDerivatives()->Print();
    event->GetAlignmentDerivatives()->Print();
    event->GetTrackPrediction()->Print();
  }
  f1.Close();
  */
}
