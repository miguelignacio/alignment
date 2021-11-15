//
//  CustomConstraints.h
//  
//
//  Created by Sebouh Paul on 11/14/21.
//
#include "TArrayD.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TEnv.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRint.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVectorD.h"
#ifndef CustomConstraints_hpp
#define CustomConstraints_hpp

#include <stdio.h>

#endif /* CustomConstraints_hpp */

/*
 * Initializes detector-specific constraints
 * to the covariance matrix
 */
void setCustomConstraints(TEnv & mEnV, TMatrixDSym & alignmentCovariance);

