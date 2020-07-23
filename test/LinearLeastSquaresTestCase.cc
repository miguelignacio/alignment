#include "LinearLeastSquaresTestCase.h"

#include "TVectorD.h"
#include "TMatrix.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"

#include "linear_least_squares.h"

#include <iostream>
#include <time.h>

using namespace std;

/** Test if two values are equal within a certain precision. */
bool equal(const double val1, const double val2, const double precision)
{
  if (val1 == 0)
    return val2 < precision;
  if (val2 == 0)
    return val1 < precision;
//   cout << "val1= " << val1 << "   val2= " << val2 << "  equality= " << TMath::Abs(2*(val1-val2)/(val1+val2)) << endl;
  return TMath::Abs(2*(val1-val2)/(val1+val2)) < precision;
}

CPPUNIT_TEST_SUITE_REGISTRATION(LinearLeastSquaresTestCase);

LinearLeastSquaresTestCase::LinearLeastSquaresTestCase()
{

}

LinearLeastSquaresTestCase::~LinearLeastSquaresTestCase()
{

}

void LinearLeastSquaresTestCase::testLinearFit()
{
  for (int dim = 4; dim < 10; dim++) {
    const int dim = 4;
    TVectorD residual(dim);
    TMatrixD covariance(dim, dim);
    TMatrixD jacobian(dim, 2);
    double x[dim];
    double y[dim];
    double erry[dim];
    
    for (int i = 0; i < dim; i++) {
      residual[i] = i + gRandom->Gaus(); // i + 1 - 2 * (i%2);
      jacobian[i][0] = 1.;
      jacobian[i][1] = i;
      covariance[i][i] = 1.;
      x[i] = i;
      y[i] = residual[i];
      erry[i] = TMath::Sqrt(covariance[i][i]);
    }
    TCanvas * c1 = new TCanvas;
    TVectorD result = linear_least_squares(residual, jacobian, covariance);
    TGraphErrors * gr = new TGraphErrors(dim, x, y, 0, erry);
    gr->Draw("Ap");
    TF1 * fit = new TF1("fit", "[0]+[1]*x", -1, dim+1);
    fit->SetParameter(0, result[0]);
    fit->SetParameter(1, result[1]);
    fit->SetLineColor(kGreen);
    fit->Draw("same");
    TF1 * user = new TF1("user", "pol1", -1, dim+1);
    gr->Fit("user", "q");
    CPPUNIT_ASSERT(equal(user->GetParameter(0), result[0], 1E-10));
    CPPUNIT_ASSERT(equal(user->GetParameter(1), result[1], 1E-10));  
    delete user;
    delete fit;
    delete gr;
    delete c1;
  }
}

void LinearLeastSquaresTestCase::testMatrixOperations()
{
  const double data[] = { 1., 2., 3., 4.};
  TMatrixD m1(2, 2, data);
  CPPUNIT_ASSERT(equal(m1[0][0], 1., 1E-10));
  CPPUNIT_ASSERT(equal(m1[0][1], 2., 1E-10));
  CPPUNIT_ASSERT(equal(m1[1][0], 3., 1E-10));
  CPPUNIT_ASSERT(equal(m1[1][1], 4., 1E-10));
  TMatrixD m2(2, 2);
  m2.Transpose(m1);
  CPPUNIT_ASSERT(m2[0][1] == m1[1][0]);
  CPPUNIT_ASSERT(m1[0][1] == m2[1][0]);
  CPPUNIT_ASSERT(m2[1][1] == m1[1][1]);
  CPPUNIT_ASSERT(m1[0][0] == m2[0][0]);
  const double cdata[] = { 1., -1, -1, 4 };
  TMatrixDSym c(2, cdata);

  TMatrixD res1 = m2*c*m1;
  CPPUNIT_ASSERT(equal(res1[0][0], 31, 1E-10));
  CPPUNIT_ASSERT(equal(res1[0][1], 40, 1E-10));
  CPPUNIT_ASSERT(equal(res1[1][0], 40, 1E-10));
  CPPUNIT_ASSERT(equal(res1[1][1], 52, 1E-10));
  TMatrixD res2 = res1;
  res1.Invert();

  TMatrixD one = res1 * res2;
  CPPUNIT_ASSERT(equal(one[0][0], 1., 1E-10));
  CPPUNIT_ASSERT(equal(one[0][1], 0., 1E-10));
  CPPUNIT_ASSERT(equal(one[1][0], 0., 1E-10));
  CPPUNIT_ASSERT(equal(one[1][1], 1., 1E-10));
}

void LinearLeastSquaresTestCase::testChi2one()
{
  const double residual[] = { 1, 0, 3, 2};
  TMatrixD covariance(4,4);
  for (int i = 0; i < 4; i++) {
    covariance[i][i] = 2.;
  }
  TMatrixD W = covariance;
  W.Invert();
  TMatrixD r(4, 1, residual);
  TMatrixD rt(1, 4, residual);
  TMatrixD chi2(1, 1); 
  chi2 = rt*W*r;
  CPPUNIT_ASSERT(equal(chi2[0][0], 7, 1E-10));
}

void LinearLeastSquaresTestCase::testChi2two()
{
  // this test fails on a statistical basis with a very low probability
  for (int i = 0; i < 100; i++) {
    gRandom->SetSeed(time(0)+1033*i);
    const int ndof = 100;
    TCanvas * c1 = new TCanvas;
    TH1F * h1 = new TH1F("h1", "Normal distribution", 100, -3, 3);
    h1->FillRandom("gaus", 10000);
    h1->Draw();
    TF1 * mygaus = new TF1("mygaus", "330./(sqrt(2))*exp(-0.5*x*x)", -2.5, 2.5);
    mygaus->Draw("same");
  
    double chi2 = 0.;
    for (int i = 0; i < ndof; i++) {
      double tmp = gRandom->Gaus(0,1);
      chi2 += tmp*tmp;
    }
    // shall be within three sigma!
    CPPUNIT_ASSERT((ndof - 4*TMath::Sqrt(2*ndof) < chi2) 
  		 && (chi2 < ndof + 4*TMath::Sqrt(2*ndof)));
    delete mygaus;
    delete h1;
    delete c1;
  }
}
