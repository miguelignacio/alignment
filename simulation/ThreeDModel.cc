// C/C++ includes
#include <iostream>

// ROOT includes
#include <TRandom.h>
#include <TMath.h>

// local includes
#include "Utilities.h"
#include "ThreeDModel.h"

using namespace std;
using namespace CLHEP;

/** \file ThreeDModel.cc
 *
 *  \brief Definition of functions and procedures for 3D manipulation with matrices.
 *
 *  Rotations are always defined to transform global coordinates to local
 *  coordinates. (GEANT convention).
 *
 */

/** Fill a rotation matrix with the Karimaki angles alpha, beta, gamma. */
void FillRotMatrixKarimaki(HepMatrix & m, double alpha, double beta, double gamma)
{
  /// rotation matrix definition from Veikko Karimaki
  HepMatrix rot_alpha(3, 3, 1);
  rot_alpha[1][1] = cos(alpha);
  rot_alpha[1][2] = -sin(alpha);
  rot_alpha[2][2] = rot_alpha[1][1];
  rot_alpha[2][1] = -rot_alpha[1][2];
  HepMatrix rot_beta(3, 3, 1);
  rot_beta[0][0] = cos(beta);
  rot_beta[0][2] = sin(beta);
  rot_beta[2][2] = rot_beta[0][0];
  rot_beta[2][0] = -rot_beta[0][2];
  HepMatrix rot_gamma(3, 3, 1);
  rot_gamma[0][0] = cos(gamma);
  rot_gamma[0][1] = -sin(gamma);
  rot_gamma[1][1] = rot_gamma[0][0];
  rot_gamma[1][0] = -rot_gamma[0][1];
//   cout << "Rotation Matrix A " << rot_alpha << endl;
//   cout << "Rotation Matrix B " << rot_beta << endl;
//   cout << "Rotation Matric C " << rot_gamma << endl;
  m = rot_gamma * rot_beta * rot_alpha;
//   cout << "Combined Rotation " << m << endl;
}

/** Extract from a rotation matrix the Karimaki angles alpha, beta, gamma. */
void GetAnglesKarimaki(const HepMatrix & m, double & alpha, double & beta, double & gamma)
{
  // compute from the given matrix m the karimaki angles
  // this only works if for all angles abs(angle) < Pi/2
//   cout << m << endl << m[2][0] << endl << m[2][1] << endl << m[1][0] << endl;
  beta = -asin(m[2][0]);
  double cbeta = cos(beta);
  alpha = asin(m[2][1]/cbeta);
  gamma = asin(m[1][0]/cbeta);
}

/** Test the Karimaki rotation matrices. */
void testAnglesKarimaki()
{
  double Pi = TMath::Pi();
  double alpha = gRandom->Uniform(Pi)-Pi/2.;
  double beta  = gRandom->Uniform(Pi)-Pi/2.;
  double gamma = gRandom->Uniform(Pi)-Pi/2.;
//   cout << "alpha = " << alpha << endl
//        << "beta  = " << beta  << endl
//        << "gamma = " << gamma << endl;
  HepMatrix m;
  FillRotMatrixKarimaki(m, alpha, beta, gamma);
  double aa, bb, cc;
  GetAnglesKarimaki(m, aa, bb, cc);
//   cout << "after FillRotMatrixKarimaki and GetAnglesKarimaki" << endl;
//   cout << "alpha = " << aa << endl
//        << "beta  = " << b << endl
//        << "gamma = " << c << endl;
  if (TMath::Abs(aa-alpha) > 1E-7) {
    ERROR("testAnglesKarimaki failed with angle alpha");
  }
  if (TMath::Abs(bb-beta) > 1E-7) {
    ERROR("testAnglesKarimaki failed with angle beta");
  }
  if (TMath::Abs(cc-gamma) > 1E-7) {
    ERROR("testAnglesKarimaki failed with angle gamma");
  }
}

/** Fill a rotation matrix with polar angles. Only two angles (theta, phi)
 * must be given, the third rotation angle is assumed to be zero.  In order to
 * get a vector in a global coordinate system with polar angle theta and
 * azimuth angle phi, one has to apply the inverse rotation to a vector along
 * z in the local coordinate frame (0, 0, 1).
 */
void FillRotMatrixPolar(HepMatrix & m, double theta, double phi)
{
  HepMatrix rot_phi(3, 3, 1);
  rot_phi[0][0] = cos(phi);
  rot_phi[0][1] = sin(phi);
  rot_phi[1][1] = rot_phi[0][0];
  rot_phi[1][0] = -rot_phi[0][1];
  HepMatrix rot_theta(3, 3, 1);
  rot_theta[0][0] = cos(theta);
  rot_theta[0][2] = -sin(theta);
  rot_theta[2][2] = rot_theta[0][0];
  rot_theta[2][0] = -rot_theta[0][2];
//   cout << "Rotation Matrix A " << rot_phi << endl;
//   cout << "Rotation Matrix B " << rot_theta << endl;
//   cout << "Rotation Matrix ATA " << rot_phi.T() * rot_phi << endl;
//   cout << "Rotation Matrix BTB " << rot_theta.T() * rot_theta << endl;
  m = rot_theta * rot_phi;
//   cout << "Combined Rotation " << m << endl;
//   cout << "Combined Rotation MTM" << m * m.T() << endl;  
}

/** Test polar angle functions (FillRotMatrixPolar) */
void testAnglesPolar()
{
  // test angles with polar coordinates
  HepMatrix m;
  double theta = gRandom->Uniform(0.2);
  double phi   = gRandom->Uniform(0.2)-0.1;
  cout << "Rotating vector to have theta = " << theta << " and phi = " << phi << endl;
  FillRotMatrixPolar(m, theta, phi);
  // vector parallel to z in local frame
  HepVector v(3); v[0] = 0; v[1] = 0; v[2] = 1;
  // transform from local coordinate system to global
  HepMatrix m1 = m.T() * v;
  HepVector z = m1;
  cout << "Vector z after rotation " << z << endl;
  cout << "z theta " << acos(z[2]) << endl;
  cout << "z phi   " << atan2(z[1], z[0]) << endl;
  HepMatrix m2 = m * m.T();
  cout << "Rotating back " << m2 << endl;
}

/** Test a matrix by generating a random matrix with polar angles and extract
    the Karimaki angles. Then the Karimaki angles are used to filled another
    matrix, and the two matrices are being compared.
*/
void testMixedAngles(int num)
{
  HepMatrix polar(3,3);
  HepMatrix karimaki(3,3);
  for (int i = 0; i < num; i++) {
    // fill polar matrix
    double theta = gRandom->Uniform(0, TMath::Pi());
    double phi   = gRandom->Uniform(0, 2*TMath::Pi());
    FillRotMatrixPolar(polar, theta, phi);
    double alpha, beta, gamma;
    // extract karimaki angles
    GetAnglesKarimaki(polar, alpha, beta, gamma);
    // fill karimaki matrix
    FillRotMatrixKarimaki(karimaki, alpha, beta, gamma);
    // compare the two
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
	double diff = TMath::Abs(polar[j][k] - karimaki[j][k]);
	if (diff > 1E-7) {
	  cout << "ERR: difference["<<j<<"]["<<k<<"] is " << diff << endl;
	  cout << "polar[][] = " << polar[j][k]
	       << " kari[][] = " << karimaki[j][k] << endl;
	}
      }
    }
  }
}

