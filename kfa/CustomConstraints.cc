//
//  CustomConstraints.cc
//  
//
//  Created by Sebouh Paul on 11/14/21.
//

#include "CustomConstraints.h"

void setCustomConstraints(TEnv & mEnv, TMatrixDSym & alignmentCovariance){
  int detector= mEnv.GetValue("customConstraintsForDetector", 0);
  
  if(detector == 1){  //CLAS12 CVT
    std::cout << "custom constraints selected for CLAS12 CVT" << std::endl;
    int nParameters=6;
    double fixationError = mEnv.GetValue("fixationError", 1E-10);
    double pairTrans = mEnv.GetValue("clas12SvtPairCovarianceTransverse", 1E-4) ;
    double pairLong = mEnv.GetValue("clas12SvtPairCovarianceLongitudinal", 1E-4) ;
    for(int i = 0; i<42; i++){
      //constraint on the space between the SVT top and bottom sensors in each module (that is, constraining that the normal component of their translation is the same).  Also, set tighter covariance of
      //first determine phi of the normal of the
      double phi = 0;
      if(i<10){ //region 1
        phi = -TMath::Pi()/2-i*2*TMath::Pi()/10.;
      } else if(i<24){ //region 2
        phi = -TMath::Pi()/2-(i-10)*2*TMath::Pi()/14.;
      } else { // region 3
        phi = -TMath::Pi()/2-(i-24)*2*TMath::Pi()/18.;
      }
      double c = TMath::Cos(phi);
      double s = TMath::Sin(phi);
      
      double deltaT2 = alignmentCovariance[i*nParameters+0][i*nParameters+0];
      
      //diagonal elements for normal and transverse:
      alignmentCovariance[i*nParameters+0][i*nParameters+0] += (c*c*fixationError + s*s*pairTrans)/4;
      alignmentCovariance[i*nParameters+1][i*nParameters+1] += (s*s*fixationError + c*c*pairTrans)/4;
      alignmentCovariance[(42+i)*nParameters+0][(42+i)*nParameters+0] += (c*c*fixationError + s*s*pairTrans)/4;
      alignmentCovariance[(42+i)*nParameters+1][(42+i)*nParameters+1] += (s*s*fixationError + c*c*pairTrans)/4;
      //now for off-diagonal elements
      alignmentCovariance[i*nParameters+0][(42+i)*nParameters+0] = deltaT2- (c*c*fixationError + s*s*pairTrans)/4;
      alignmentCovariance[(42+i)*nParameters+0][i*nParameters+0] = deltaT2- (c*c*fixationError + s*s*pairTrans)/4;
      alignmentCovariance[i*nParameters+1][(42+i)*nParameters+1] = deltaT2- (s*s*fixationError + c*c*pairTrans)/4;
      alignmentCovariance[(42+i)*nParameters+1][i*nParameters+1] = deltaT2- (s*s*fixationError + c*c*pairTrans)/4;
      //mix Tx and Ty
      alignmentCovariance[i*nParameters+0][i*nParameters+1] = c*s*(fixationError-pairTrans)/4;
      alignmentCovariance[i*nParameters+1][i*nParameters+0] = c*s*(fixationError-pairTrans)/4;
      alignmentCovariance[(i+42)*nParameters+0][(i+42)*nParameters+1] = c*s*(fixationError-pairTrans)/4;
      alignmentCovariance[(i+42)*nParameters+1][(i+42)*nParameters+0] = c*s*(fixationError-pairTrans)/4;
      //mix Tx and Ty of opposite sensors
      alignmentCovariance[i*nParameters+0][(i+42)*nParameters+1] = -c*s*(fixationError-pairTrans)/4;
      alignmentCovariance[(i+42)*nParameters+1][i*nParameters+0] = -c*s*(fixationError-pairTrans)/4;
      alignmentCovariance[(i+42)*nParameters+0][i*nParameters+1] = -c*s*(fixationError-pairTrans)/4;
      alignmentCovariance[i*nParameters+1][(i+42)*nParameters+0] = -c*s*(fixationError-pairTrans)/4;
      
      
      // z constraint
      alignmentCovariance[i*nParameters+2][i*nParameters+2] += (pairLong)/4;
      alignmentCovariance[(42+i)*nParameters+2][(42+i)*nParameters+2] += (pairLong)/4;
      alignmentCovariance[i*nParameters+2][(42+i)*nParameters+2] = deltaT2-pairLong/4;
      alignmentCovariance[(42+i)*nParameters+2][i*nParameters+2] = deltaT2-pairLong/4;
      
      
      
      //constraint that the SVT top and bottom sensors in each module
      //have the same rotation
      for(int j = 3; j<6; j++){
        alignmentCovariance[i*nParameters + j][(42+i)*nParameters + j] = (alignmentCovariance[(42+i)*nParameters + j][(42+i)*nParameters + j]+alignmentCovariance[i*nParameters + j][i*nParameters + j]-fixationError)/2;
        alignmentCovariance[(42+i)*nParameters + j][i*nParameters + j]=alignmentCovariance[i*nParameters + j][(42+i)*nParameters + j];
      }
    }
  }
  
}
