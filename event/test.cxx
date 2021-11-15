#include <iostream>

#include "AlignEvent.h"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

using namespace std;

int main(int argc, char * argv[]) 
{
  const int nEvents = 10;
   // initialize ROOT
  TROOT root("root", "root");
  std::cout << "hello world" <<std::endl;
  
  // create event
  AlignEvent * event = new AlignEvent;
  
  std::cout << "hello world" <<std::endl;
  
}
