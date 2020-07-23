#include "linear_least_squares.h"

#include <iostream>

using namespace std;

#include "TRandom.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"

TVectorD linear_least_squares(const TVectorD & residual, const TMatrixD & Jacobian, const TMatrixD & covariance)
{
  // compute weight matrix
  double determinant;
  const int npar = Jacobian.GetNcols();
  TMatrixD parameter_covariance(npar, npar);
  TMatrixD W = covariance;
  TVectorD result(npar);
  W.Invert(& determinant);
  if (determinant == 0) {
    cout << "ERR: Determinant is zero / cannot solve problem" << endl;
    return result;
  }
  // cout << "Weight matrix"; W.Print();

  // compute chi2 before fit
  TMatrixD r(residual.GetNrows(), 1, residual.GetMatrixArray());
  TMatrixD rt(1, residual.GetNrows(), residual.GetMatrixArray());
  // cout << "residual = "; residual.Print();
  // cout << "r = "; r.Print();
  // cout << "rt = "; rt.Print();
  TMatrixD chi2(1,1);
  chi2 = rt*W*r;
  // cout << "chi2 before fit: " << chi2[0][0] << endl;

  // perform fit, compute fit parameters
  TMatrixD JT(npar, Jacobian.GetNrows());
  JT.Transpose(Jacobian);
  parameter_covariance = JT*W*Jacobian;
  // cout << "parameter_covariance matrix"; parameter_covariance.Print();
  parameter_covariance.Invert(& determinant);
  if (determinant == 0) {
    cout << "ERR: Determinant is zero / cannot solve problem" << endl;
    return result;
  }
  // cout << "parameter covariance"; parameter_covariance.Print();
  // cout << "result init:"; result.Print();
  result = parameter_covariance*JT*W*residual;
  // cout << "result comp:"; result.Print();

  // compute corrected residual
  TVectorD residual_new = residual - Jacobian*result;

  // compute chi2 after fit
  TMatrixD r2(residual_new.GetNrows(), 1, residual_new.GetMatrixArray());
  TMatrixD rt2(1, residual_new.GetNrows(), residual_new.GetMatrixArray());
  // cout << "residual = "; residual_new.Print();
  // cout << "r2 = "; r2.Print();
  // cout << "rt2 = "; rt2.Print();
  chi2 = rt2*W*r2;
  // cout << "chi2 (after fit) = " << chi2[0][0]
  //      << " ndof = " << residual.GetNrows()-npar << endl;

  // compute covariance matrix of parameters
  return result;
}
