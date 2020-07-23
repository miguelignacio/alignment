/**
 *
 * \mainpage Kalman Filter Alignment
 *
 * This is a toy model with which you can simulate Kalman Filter Alignment. It
 * consists of a simulation package, the kalman alignment algorithm, and a
 * validation package. Multiple scattering is simulated, and it is possible to
 * generate particles according to the cosmic muon spectrum.
 *
 * Build instructions:
 *
 * This package relies on two external libraries, namely CLHEP
 * (http://proj-clhep.web.cern.ch/proj-clhep/) and ROOT
 * (http://root.cern.ch/). You will first need to install these packages.
 * Adjust the SConstruct and the event/Makefile and util/Makefile accordingly.
 *
 * Currently two builder programs are used: SCons and make. First, you need to
 * build two shared libraries. Go to the "event/" subdirectory and issue
 * "make". Add the "event/" directory to your LD_LIBRARY_PATH. Go to the
 * "util/" subdirectory and issue "make". Add the "util/" directory to your
 * LD_LIBRARY_PATH. Go back to the main directory (where the Doxyfile and the
 * SConstruct file reside). Run the SCons program by typing "scons" at the
 * prompt. This will build all other programs. Run the simulation in the main
 * directory with the command
 * 
 *
 * simulation/simulation cfg/simulation.cfg
 *
 * The program will simulate 10000 tracks in a misaligned detector,
 * reconstruct the tracks, and show a display of all tracks. A file with all
 * information necessary for the alignment of the detector is created, by
 * default called "alignment_data.root". Open ROOT, load the shared library
 * and open the alignment_data.root file. Then you can browse through it. An
 * alternative is to use the event/event program to display the first 10 data
 * records.
 *
 * After running the simulation, it is time for alignment. Run the alignment
 * program with
 *
 * kfa/align cfg/align.cfg
 *
 * This runs the Kalman Alignment algorithm over the previously created data
 * file. After finishing, it has created another file named
 * "alignment_result.root". This file contains for each alignment parameter
 * its evolution, i.e. its value for each processed track. These histograms
 * are named "hAliParXXX" where XXX is the number of the alignment parameter.
 * In addition, the final values of the alignment parameters are saved in the
 * histogram hLimits. The values in each bin correspond to the alignment
 * parameter after processing all tracks. The bin number corresponds to the
 * alignment parameter number.
 *
 * Finally, you can validate whether the alignment procedure worked. Run the
 * program
 *
 * validation/validation cfg/validation.cfg
 * 
 * This will start the data-driven validation. The validation program will
 * reprocess all tracks with the previously determined alignment parameters
 * and compare the track chi^2 value before and after the alignment. The
 * program creates three plots: A raw plot named "validation_1.pdf" and more
 * refined plots "validation_2.pdf" and "validation_3.pdf". Looking at the
 * latter ones is recommended. They shows a comparison between the track chi^2
 * before alignment and after alignment. The file "validation_2.pdf" shows the
 * track chi^2 directly. After alignment, the mean chi^2 should correspond to
 * the number of degrees of freedom of each track. For each crossed layer, the
 * track has two hits. E.g. for the default detector with six layers there are
 * twelve hits. Since the track has four parameters (no magnetic field is
 * simulated), the number of degrees of freedom for the track is 2*6-4 = 8. So
 * you should expect that after alignment the mean chi^2 is approximately
 * eight. In practive it will be smaller since some tracks do not cross some
 * detectors, so the number of hits is on average smaller than 2*6. You can
 * see the number of deegrees of freedom in the "validation_1.pdf" histogram.
 * The third histogram, "validation_3.pdf", shows the so-called "normalized"
 * chi^2-distribution, i.e. the value "chi^2/ndof" where ndof is the number of
 * degrees of freedom as above. After alignment it is expected to be close to
 * unity.
 * 
 *
 * You can configure the simulation, alignment and validation packages by
 * editing the corresponding files in the cfg/ subdirectory.
 *
 * Have fun!
 *
 * 
 * Copyright &copy; 2005-2011 Martin Weber, Daniel Sprenger
 *
 */

/** \file simulation.cxx 
 *
 * \brief This is the main file for simulation and track reconstruction.
 *
 */

// C/C++ includes
#include <iostream>
#include <iomanip>
#include <string>

// STL includes
#include <vector>

// ROOT includes
#include "TROOT.h"
#include "TRint.h"
#include "TRandom.h"
#include "TBenchmark.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TEnv.h"
#include "TMath.h"
#include "TStyle.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TFile.h"
#include "TTree.h"
#include "TArray.h"

// CLHEP includes
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"

// local includes
#include "Utilities.h"
#include "Constants.h"
#include "Det.h"
#include "Detector.h"
#include "Track.h"
#include "ThreeDModel.h"
#include "Hit.h"
#include "KalmanFilterTracking.h"
#include "Visualization.h"
#include "AlignEvent.h"
#include "KalmanFilterAlignmentInputProvider.h"

using namespace std;
using namespace CLHEP;

//////////////////////////////////////////////////////////////////////
// Global variables

/// Global Visualization object
Visualization gVis;

/// Global configuration for minimum energy (e.g. lead shield simulation)
double gMinimumEnergy;

/// How to simulate detector resolution
enum DetectorResolutionFunction { kDelta, /// perfect detector resolution
				  kUniform, /// simulate uniform resolution
				  kGaus }; /// simulate gaussian resolution

/** \brief Generate a track, distributed according to the cosmic muon spectrum
 *  over the whole detector area.
 *
 *  The muon momentum distribution is taken from
 *  http://arxiv.org/pdf/hep-ph/0103322.pdf and has been fitted by a third
 *  order polynomial. Generation is only done in a cone with an opening angle
 *  of 20 degrees, according to the Caprice measurements.
 *
 */
Track generate_cosmic_track()
{
  static TF1 * fp = 0;
  static TF1 * ftheta = 0;
  assert(gDetector != 0);
  const double width  = gDetector->GetWidth();
  const double length = gDetector->GetLength();
  const double height = gDetector->GetHeight();
  if (fp == 0) {
    assert((gMinimumEnergy > 0.) && gMinimumEnergy <= 100.);
    fp = new TF1("cosmic_muon_momentum", "(x^(-1.5))*10^(1.189+0.9604*log10(x)-1.114*log10(x)^2+0.1775*log10(x)^3)", gMinimumEnergy, 100);
//     const double layer_distance = gDetector->GetLayerDistance();
//     // a track should at least cross two layers
//     const double theta_max = TMath::ATan2(TMath::Max(width, length), layer_distance);
    // a track should at least have the possibility to hit all layers
    const double theta_max_acceptance = TMath::ATan2(TMath::Max(width, length), height);
    INFO("Max track angle acceptance " << theta_max_acceptance << " rad" );
    // model validity for only 10 degrees...
    double theta_max = 10./180.*TMath::Pi();
    if (theta_max < theta_max_acceptance) {
      WARNING("Detector acceptance is larger than cosmic ray generator allows (" << theta_max << " rad)");
      WARNING("Generating tracks only within valid CAPRICE measurement range");
//       WARNNIG("I will continue with the given detector acceptance " << theta_max_acceptance << " rad");
//       theta_max = theta_max_acceptance;
    }
    else {
      theta_max = theta_max_acceptance;
    }
    // function for angular spectrum
    ftheta = new TF1("cosmic_muon_angle", "cos(x)*cos(x)", 0, theta_max);
//     new TCanvas("cf", "cosmic function");
//     ftheta->Draw();
//     fp->Draw();
  }
  double momentum = fp->GetRandom();
  double theta    = ftheta->GetRandom();
  double phi      = gRandom->Uniform(-TMath::Pi(), TMath::Pi());

  LOG(5, "track momentum " << momentum 
      << " theta " << theta << " phi " << phi);
  HepVector origin;
  setXYZ(origin, 
	 gRandom->Uniform(width), 
	 gRandom->Uniform(length),
	 gRandom->Uniform(height));
  HepVector destination;
  setXYZ(destination,
	 origin[0]+TMath::Cos(phi)*TMath::Sin(theta)*height,
	 origin[1]+TMath::Sin(phi)*TMath::Sin(theta)*height,
	 height);
  Track t(origin, destination, momentum);
  return t;
}

/** \brief Generate a track, uniformely distributed over the whole detector area.
 *
 *  \param momentum is the track momentum. It default value is 4 GeV, which is
 *  the mean energy of the cosmic muon momentum distribution.
 *
 */
Track generate_uniform_track(double momentum = 4)
{
  assert(gDetector != 0);
  double width = gDetector->GetWidth();
  double length = gDetector->GetLength();
  double height = gDetector->GetHeight();
  HepVector origin;
  setXYZ(origin, 
	 gRandom->Uniform(width), 
	 gRandom->Uniform(length),
	 0.);
  HepVector destination;
  setXYZ(destination, 
	 gRandom->Uniform(width),
	 gRandom->Uniform(length),
	 height);
  Track t(origin, destination, momentum);
  return t;
}


/** \brief Check the measurement model
 *
 * This function compares the predicted hit (intersection of track with the
 * nominal detector position) with the calculation from the measurement model.
 * If the true track parameters and the true misalignment is known, the
 * measurement model allows to compute the true position from the predicted
 * position.
 *
 * The values reported should be identical. However, this is only true if
 * there has been no simulation of detector resolution and no multiple
 * scattering effects are being computed. If those are incorporated,
 * differences are expected.
 *
 */
void check_measurement_model(const Track & track, const vector<Hit *> hitvector)
{
  for (vector<Hit *>::const_iterator it = hitvector.begin();
       it != hitvector.end(); it++) {
    const Det & adet = (*it)->GetDet();
    // get track state in nominal position
    HepVector * r = track.IntersectionWithPlane(adet.GetNominal());
    if (r == 0) {
      WARNING("hit in misaligned detector but not in nominal detector");
      WARNING("     (this can happen if the track is at the very edge of the detector)");
      continue;
    }
    // let's do it the fast way... (least number of matrix multiplications)
    // hit in global coordinates
    cout << "r" << *r << endl;
    // translate hit to local frame
    HepVector p = adet.GetNominal().ToLocal(*r);
    cout << "p" << p << endl;
    // translate hit in misaligned frame
    HepVector qprime = adet.GetMisalign().ToLocal(*r);
    cout << "qprime = " << qprime << endl;
    // now project along trajectory. Use misaligned frame.
    HepVector t = track.GetLocalParameters(adet.GetNominal());
    HepVector n(3, 1);
    n[0] = t[0];
    n[1] = t[1];
    cout << "n = " << n << endl;
    // rotate in correct frame
    n = adet.GetMisalignDelta().GetRotation() * n;
    // adjust length
    n = n * (1./n[2]);
    cout << "nprime = " << n << endl;
    cout <<  qprime[2] << endl;
    HepVector m = qprime - qprime[2] * n;
    cout << "m = " << m << endl;
    cout << "should be equal to " << (*it)->GetPosition();
  }
}

/** \brief Calculate the projected mean angular deflection due to multiple
 *  scattering of a particle in material.
 *
 *  Projected means that the angle is not the angle in space, but rather a
 *  one-dimensional projection on one axis of a plane that is perpendicular to
 *  the particle direction. To obtaint the spatial angle, one has to use twice
 *  the returned value.
 *
 *  \param beta Velocity/c
 *  \param z    Charge of particle)
 *  \param p    Momentum, must be given in MeV
 *  \param x    Matter thickness in cm
 *  \param xi   Radiation length of material in cm
 * 
 */
double ms_angle(double beta, int z, double p, double x, double xi)
{
  return 13.6 / (beta*p) * z * TMath::Sqrt(x/xi)*(1+0.038*TMath::Log(x/xi));
}


/** \brief Simulation of multiple scattering. 
 *
 * The track is scattered at the given detector. The offset through scattering
 * in the material is neglected.
 *
 * \param t    The track before (and after) multiple scattering. The track
 *             parameters are modified according to the new direction.
 *
 * \param adet Detector in which the multiple scattering occurs.
 *
 * \param hit  Must be the hit before applying any other smearing, because it
 *             is taken as the starting point for the resulting track!
 *
 * The routine estimates the radiation length according to the path length of
 * the particle track through the silicon detector. 
 */
void simulate_multiple_scattering(Track & t, const Det & adet, const Hit * hit)
{
  double p = t.GetMomentum();
  double m = t.GetMass();
  int    c = TMath::Abs(t.GetCharge());
  double E = TMath::Sqrt(p*p + m*m);
  LOG(5, "E = " << E << " GeV");
  double gamma = E/m;
  LOG(5, "gamma = " << gamma);
  double beta = TMath::Sqrt(1-1/(gamma*gamma));
  LOG(5, "beta = " << beta);

  HepVector origin = hit->GlobalTruePosition();
  HepVector destination = t.GetDestination();
  LOG(5, "Original destination: " << destination);
  HepVector direction = destination - origin;
  double distance = direction.norm();

  // Calculate the approximate distance the particle travels through the detector.
  HepVector normal;
  setXYZ(normal, 0, 0, 1);
  double ct = dot(direction, adet.GetMisalign().GetRotation().T()*normal)/distance;
  LOG(5, "ct = " << ct);

  // Compute multiple scattering deflection angle
  double angle = 2.*ms_angle(beta, 
			     c, 
			     p*1000,
			     adet.GetThickness()/10./ct,
			     adet.GetRadiationLength());
  LOG(5, "angle = " << angle);
  // transform angle into an offset at the destination plane. This is not the
  // multiple scattering offset.
  double offset   = TMath::Abs(gRandom->Gaus(0, distance*TMath::Tan(angle)));
  LOG(5, "track offset = " << offset);
  // create two arbitray perpendicular vectors to direction to implement
  // deflection
  double x = direction[0]; 
  double y = direction[1];
  double z = direction[2];
  HepVector dir1;
  setXYZ(dir1, 0, z, -y);
  HepVector dir2; 
  setXYZ(dir2, -y*y-z*z, x*y, x*z);
  // change vectors to unit size
  dir1 /= dir1.norm();
  dir2 /= dir2.norm();
  // take flat random combination of the two vectors
  double phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
  destination += offset*(TMath::Sin(phi)*dir1 + TMath::Cos(phi)*dir2);
  LOG(5, "MS destination: " << destination);
  Track scattered_track(origin, destination, p);
  t = scattered_track;
}

/** \brief Simulation step: Simulate hits in real (misaligned) detector.
 *									
 * Compute hits in all detectors for the given track. Tracks are smeared with
 * the detector resolution if asked for. */
void simulate_hits(Track & t, 
		   DetectorResolutionFunction smear = kUniform, 
		   bool multiple_scattering = false)
{
  // get detector layout
  assert(gDetector != 0);
  int nlayer  = gDetector->GetNLayer();
  int nrod    = gDetector->GetNRod();
  int nmodule = gDetector->GetNModule();

  // loop over layers in order to calculate hits
  for (int l = 0; l < nlayer; l++) {
    bool hit_found = false;
    for (int r = 0; r < nrod; r++) {
      if (hit_found)
	break;
      for (int m = 0; m < nmodule; m++) {
	const Det & adet = gDetector->GetDet(l, r, m);
	// we have to calculate the hit in the misaligned frame since the
        // detector is there in reality and we are computing real hits here,
        // not predicted ones!
	Hit * hit = t.MisalignedIntersection(adet);
	if (!hit)
	  continue;
	// cout << "Hit true position: " << hit->GlobalTruePosition();
	// simulate multiple scattering with given distance
	if (multiple_scattering) {
	  simulate_multiple_scattering(t, adet, hit);
	}
	if (smear != kDelta) {
	  // smear hit with detector resolution
	  double pitch_x = adet.GetPitchX();
	  double pitch_y = adet.GetPitchY();
	  // we need to make sure that the smeared hit is inside the detector
	  // bounds. However, this introduces a small bias because simulated
	  // smeared hits are on average closer to the true position. Ideally
	  // one should assign different errors to hits at the detector edge to
	  // avoid this problem.
	  while (true) {
	    HepVector position = hit->GetPosition();
	    if (smear == kUniform) {
	      // using "+=" in the next line leads to a random walk if the hit
	      // is outside bounds, therefore the while loop can take a long
	      // time to finish. One could improve here.
	      position[0] += gRandom->Uniform(-pitch_x/2., pitch_x/2.);
	      position[1] += gRandom->Uniform(-pitch_y/2., pitch_y/2.);
	    } else if (smear == kGaus) {
	      // using "+=" in the next line leads to a random walk if the hit
	      // is outside bounds, therefore the while loop can take a long
	      // time to finish. One could improve here.
	      position[0] += gRandom->Gaus(0, pitch_x/TMath::Sqrt(12.)); 
	      position[1] += gRandom->Gaus(0, pitch_y/TMath::Sqrt(12.)); 
	    } else
	      THROW(std::string("Unknown detector resolution function"));
	    // make sure the hit is inside bounds
	    if (hit->GetDet().InsideBounds(position)) {
	      // ok, is inside, set position and break loop
	      hit->SetPosition(position);
	      break;
	    }
	  }
	}
	t.GetHits().push_back(hit);
	// cout << "Hit global measured position: " << hit->GlobalMeasuredPosition();
	// break loop, since we found a hit
	hit_found = true;
	break;
      }
    }
  }
  // cout << "Found " << t.GetHits().size() << " hits on track." << endl;
}

/** \brief Test transformation between global and local track parameters.
 */
void check_track_parameter_transformation(Track t)
{
  assert(gDetector != 0);
  int nlayer = gDetector->GetNLayer();
  
  for (int i = 0; i < nlayer; i++) {
    const Det & adet = gDetector->GetDet(i, 0, 0);
    cout << "Original track: " << t.GetOrigin() << t.GetDestination() << endl;
    HepVector v = t.GetGlobalParameters();
    cout << "Global parameters " << v << endl;
    v = t.GetLocalParameters(adet.GetMisalign());
    cout << "Local parameters in det["<<i<<"][0][0]" << v << endl;
    t.SetLocalParameters(adet.GetMisalign(), v);
    cout << "After setting those parameters"
	 << t.GetOrigin() << t.GetDestination() << endl;
  }
}

/** \brief Test transformation between global and local track covariance
 */
void check_track_covariance_transformation()
{
  assert(gDetector != 0);
  int nlayer  = gDetector->GetNLayer();

  for (int i = 0; i < 100; i++) {
    HepSymMatrix covariance(4,0);
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
	covariance[j][k] = gRandom->Uniform(-1., 1.);
      }
    }
    const Det & adet = gDetector->GetDet(i%nlayer, 0, 0);
    Track t; 
    t.SetGlobalCovariance(covariance);
    HepSymMatrix local = t.GetLocalCovariance(adet.GetMisalign());
    t.SetLocalCovariance(adet.GetMisalign(), local);
    HepSymMatrix back = t.GetGlobalCovariance();
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
	double diff = covariance[j][k]-back[j][k];
	if (TMath::Abs(diff)>1E-7) {
	  cout << "i=" << i << " j=" << j << " k=" << k << " diff= " << diff << endl;
	}
      }
    }
  }
}

void check_track_linearization()
{
  assert(gDetector != 0);
  int nlayer  = gDetector->GetNLayer();
  int nrod    = gDetector->GetNRod();
  int nmodule = gDetector->GetNModule();

  // the track we start with
  Track t = generate_uniform_track();
  const int nPoints = 100; ///< plot nPoints points
  const double range = 1;  ///< plot in [-range,range) around current track parameters
  const double delta = 2*range / nPoints;

  /// create canvas for plotting
  TCanvas * ctl = new TCanvas("ctl", "Track linearization", 10, 10, 1024, (int) (1024./TMath::Sqrt(2.)));
  ctl->Divide(4,4);
  ctl->Draw();

  /// This is where the track linearization takes place
  HepVector pglobal0 = t.GetGlobalParameters();
  cout << "global parameters" << pglobal0;

  /// loop over all detectors
  for (int l = 0; l < nlayer; l++) {
    for (int r = 0; r < nrod; r++) {
      for (int m = 0; m < nmodule; m++) {

	/// get detector
	const Det & adet = gDetector->GetDet(l, r, m);

	/// compute track parameters at detector location
	const ReferenceFrame & frame = adet.GetMisalign();

	// reset track for next detector
	t.SetGlobalParameters(pglobal0); 
	HepVector plocal0 = t.GetLocalParameters(frame);

	/// loop over all four parameters and fill extrapolation histograms
	for (int i = 0; i < 4; i++) {
	  TH1D * h1[4];
	  TH1D * h2[4];
	  /// create histograms
	  for (int k = 0; k < 4; k++) {
	    h1[k] = new TH1D(Form("hlocal_%d_%d", i, k), Form("Parameter (linearized)%d %d", i, k), nPoints, -range, range);
	    h2[k] = new TH1D(Form("hlocal0_%d_%d", i, k), Form("Parameter (exact)  %d %d", i, k), nPoints, -range, range);
	  }
	  for (int j = 0; j < nPoints; j++) {
	    HepVector pglobal = pglobal0;
	    double x = delta*j-range+(delta/100.); ///< add small step to avoid strange binning effect (center of bin)
	    pglobal[i] += x;
	    HepMatrix Jacobian = t.Jacobian(frame);
	    t.SetGlobalParameters(pglobal);
	    HepVector extrapolation = plocal0 + Jacobian * (pglobal - pglobal0);
	    HepVector plocal = t.GetLocalParameters(frame); 
// 	    if ((i == 3) && (TMath::Abs(x) < 0.001)) {
// 	      cout << "x = " << x << endl;
// 	      cout << "plocal0 " << plocal0 << "plocal" << plocal;
// 	      cout << "pglobal0 " << pglobal0 << "pglobal" << t.GetGlobalParameters();
// 	      t.SetLocalParameters(frame, plocal0);
// 	      cout << "back side " << t.GetGlobalParameters();
// 	    }
	    /// fill histograms
	    for (int k = 0; k < 4; k++) {
	      h1[k]->Fill(x, extrapolation[k]-plocal0[k]);
	      h2[k]->Fill(x, plocal[k]-plocal0[k]);
	    }
	  }
	  // draw histograms
	  for (int k = 0; k < 4; k++) {
	    ctl->cd(i*4+k+1);
	    h1[k]->SetLineColor(kBlue);
	    h1[k]->Draw();
	    h2[k]->SetLineColor(kRed);
	    h2[k]->Draw("same");
	  }
	}
	ctl->Update();
	char buffer;
	cin >> buffer;
      }
    }
  }
}

/** \brief Compute tracking pulls, difference to truth dividided by tracking errors.
 *
 * This function uses the known truth, i.e. the true track parameters, and
 * compares the reconstructed ones to them, dividing by the track errors. This
 * "pull" distribution should be, in case of correct parameters and error
 * estimates, have mean zero and RMS one. With gaussian input errors
 * (i.e. gaussian hit errors), it should follow a gaussian
 * distribution. However, hits are distributed according to uniform
 * distribution, therefore the shape is not gaussian (but mean should still be
 * zero and RMS one).
 */
void fill_tracking_pulls(const Track & reco, const Track & truth)
{
  static TH1F * hp1 = 0;
  static TH1F * hp2 = 0;
  static TH1F * hp3 = 0;
  static TH1F * hp4 = 0;
  if (hp1 == 0 || hp2 == 0 || hp3 == 0 || hp4 == 0) {
    hp1 = new TH1F("hp1", "Pull of track parameter tx", 100, -3., 3.);
    hp2 = new TH1F("hp2", "Pull of track parameter ty", 100, -3., 3.);
    hp3 = new TH1F("hp3", "Pull of track parameter x", 100, -3., 3.);
    hp4 = new TH1F("hp4", "Pull of track parameter y", 100, -3., 3.);
  }
  HepVector diff = reco.GetGlobalParameters()-truth.GetGlobalParameters();
  hp1->Fill(diff[0]/TMath::Sqrt(reco.GetGlobalCovariance()[0][0]));
  hp2->Fill(diff[1]/TMath::Sqrt(reco.GetGlobalCovariance()[1][1]));
  hp3->Fill(diff[2]/TMath::Sqrt(reco.GetGlobalCovariance()[2][2]));
  hp4->Fill(diff[3]/TMath::Sqrt(reco.GetGlobalCovariance()[3][3]));
}

/** \brief Make a histogram of track chi2, chi2/ndof and chi2-Probability. 
 */
void fill_track_chi2(const Track & reco)
{
  static TH1F * hchi2 = 0;
  static TH1F * hchi2ndof = 0;
  static TH1F * hchi2prob = 0;
  // if either histograms or canvas is not there, create them
  if (hchi2 == 0 || hchi2ndof == 0 || hchi2prob == 0) {
    hchi2 = new TH1F("hchi2", "track #chi^{2}", 15*gDetector->GetNLayer(), 0, 5*gDetector->GetNLayer());
    hchi2ndof = new TH1F("hchi2ndof", "track #chi^{2}/ndof", 100, 0, 10);
    hchi2prob = new TH1F("hchi2prob", "track #chi^{2} probability", 100, 0, 1);
  }
  hchi2->Fill(reco.GetChi2());
  hchi2ndof->Fill(reco.GetChi2()/reco.GetNdof());
  hchi2prob->Fill(TMath::Prob(reco.GetChi2(), reco.GetNdof()));
}

/** Test HepMatrix.sub() */
void subtest()
{
  HepMatrix A(36,36);
  for (int i = 0; i < 36; i++) {
    for (int j = 0; j < 36; j++) {
      A[i][j] = 36*i+j;
    }
  }
  cout << A << A.sub(1, 6, 1, 6) << std::endl;
}

int main(int argc, char * argv[])
{
  cout << 
    "Kalman Alignment - track based alignment of high-energy particle physics detectors.\n"
    "(simulation package)\n"
    "\n"
    "Copyright (C) 2005-2011 Martin Weber\n"
    "\n"
    "This program is free software: you can redistribute it and/or modify\n"
    "it under the terms of the GNU General Public License as published by\n"
    "the Free Software Foundation, either version 3 of the License, or\n"
    "(at your option) any later version.\n"
    "\n"
    "This program is distributed in the hope that it will be useful,\n"
    "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
    "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
    "GNU General Public License for more details.\n"
    "\n"
    "You should have received a copy of the GNU General Public License\n"
    "along with this program.  If not, see <http://www.gnu.org/licenses/>.\n";
  TROOT ROOT("simulation", "simulation");

  // Create interactive interface
  int argROOT = 1;
  TRint *theApp = new TRint("Simulation for Kalman Filter Alignment", &argROOT, argv, NULL, 0);

  // check command line arguments
  if (argc != 2) {
    cout << "Usage:" << endl
	 << endl
	 << "simulation ConfigFile" << endl
	 << endl
	 << "ConfigFile: Contains the configuration for detector layout, misalignment," << endl
	 << "            simulation, verbosity and more." << endl
	 << endl;
    exit(1);
  }

  // setup draw options
  gStyle->SetOptStat(1111111111);

  // read config file
  TEnv mEnv(argv[1]);

  // configure values from config file
  const char * DetectorConfiguration  = mEnv.GetValue("DetectorConfiguration", "CosmicRack.cfg");
  const char * AlignmentDataFileName  = mEnv.GetValue("AlignmentDataFileName", "alignment_data.root");
  double       PositionAlignmentError = mEnv.GetValue("PositionAlignmentError", 1.);
  double       RotationAlignmentError = mEnv.GetValue("RotationAlignmentError", 0.050);
  int          NumberOfTracks         = mEnv.GetValue("NumberOfTracks", 1);
  int          MinimumNumberOfHits    = mEnv.GetValue("MinimumNumberOfHits", 2);
  const char * ParticleGenerator      = mEnv.GetValue("ParticleGenerator", "Cosmic");
  double       GunMomentum            = mEnv.GetValue("Gun.Momentum", 4.);
  const char * DetectorResolution     = mEnv.GetValue("DetectorResolution", "Uniform");
  bool         MultipleScattering     = mEnv.GetValue("MultipleScattering", true);
               gMinimumEnergy         = mEnv.GetValue("Cosmic.MinimumEnergy", 0.1);
  gLogLevel                           = mEnv.GetValue("LogLevel", 3);

  INFO("DetectorConfiguration  = " << DetectorConfiguration);
  INFO("PositionAlignmentError = " << PositionAlignmentError);
  INFO("RotationAlignmentError = " << RotationAlignmentError);
  INFO("NumberOfTracks         = " << NumberOfTracks);
  INFO("MinimumNumberOfHits    = " << MinimumNumberOfHits);
  INFO("ParticleGenerator      = " << ParticleGenerator);
  INFO("Gun.Momentum           = " << GunMomentum);
  INFO("DetectorResolution     = " << DetectorResolution);
  INFO("MultipleScattering     = " << MultipleScattering);
  INFO("Cosmic.MinimumEnergy   = " << gMinimumEnergy);
  INFO("LogLevel               = " << gLogLevel);

  // check parameters
  DetectorResolutionFunction smear = kUniform;
  if (!strcmp(DetectorResolution, "Off") || 
      !strcmp(DetectorResolution, "off") || 
      !strcmp(DetectorResolution, "OFF") || 
      !strcmp(DetectorResolution, "false") ||
      !strcmp(DetectorResolution, "False") || 
      !strcmp(DetectorResolution, "FALSE") || 
      !strcmp(DetectorResolution, "0")) {
    smear = kDelta;
  } else if (!strcmp(DetectorResolution, "Uniform") || 
	     !strcmp(DetectorResolution, "UNIFORM") || 
	     !strcmp(DetectorResolution, "uniform")) {
    smear = kUniform;
  } else if (!strcmp(DetectorResolution, "Gaus") || 
	     !strcmp(DetectorResolution, "GAUS") || 
	     !strcmp(DetectorResolution, "gaus")) {
    smear = kGaus;
  } else {
    ERROR("Wrong value \"" << DetectorResolution <<
	  "\" given for parameter \"DetectorResolution\""
	  "\nValid entries are \"Off\", \"Uniform\" or \"Gaus\".");
    return 1;
  }

  gDetector = new Detector(DetectorConfiguration);

  // subtest(); // there was a bug once in CLHEP such that it did not work

  // for (int i = 0; i < 10000; i++) {
  //   testAnglesKarimaki();
  // }

  // check_track_linearization();

  gBenchmark->Start("simulation");

  int nlayer  = gDetector->GetNLayer();
  int nrod    = gDetector->GetNRod();
  int nmodule = gDetector->GetNModule();

  // Add detector to visualization
  gVis.AddDetector(*gDetector);

  // create tracking class
  KalmanFilterTracking tracking;

  //////////////////////////////////////////////////////////////////////
  // Data I/O

  // Open file for storing information
  TFile * alignment_data = new TFile(AlignmentDataFileName, "RECREATE");
  if (alignment_data == 0 || alignment_data->IsOpen() != kTRUE) {
    ERROR("Could not open alignment data file " << AlignmentDataFileName);
  }
  // This simulation has nlayer*nrod*nmodule alignables. 
  // We align all six degrees of freedom.
  AlignInfo * myInfo = new AlignInfo(nlayer*nrod*nmodule, 6);
  myInfo->Write();
  // create alignment "event" which records all information from one track
  AlignEvent * myEvent = new AlignEvent;
  // we assume always run number one
  myEvent->SetRun(1);
  // create a tree where the alignment info is being stored
  TTree * AlignTree = new TTree("AlignTree", "Alignment data", 0);
  AlignTree->Branch("AlignEvent", & myEvent, 64000, 0);
  // auto save every MB
  AlignTree->SetAutoSave(1000000);

  //////////////////////////////////////////////////////////////////////
  // Simulation loop
  int i = 0;
  while (i < NumberOfTracks) {
    for (int j = 0; j < 44; j++)
      cout << "\b";
    cout << setw(5) << i << flush;

    // generate a track
    Track simTrack;
    if (!strcmp(ParticleGenerator, "Cosmic")) {
      simTrack = generate_cosmic_track();
    }
    else if (!strcmp(ParticleGenerator, "Gun")) {
      simTrack = generate_uniform_track(GunMomentum);
    }
    else {
      THROW(std::string("Do not know how to generate a \"")+std::string(ParticleGenerator)+std::string("\" track!"));
    }
    // check_track_parameter_transformation(t);
    // check_track_covariance_transformation();
    // simulate the detector with detector resolution and multiple scattering effects
    // cout << "True track parameters: " << simTrack.GetGlobalParameters();
    simulate_hits(simTrack, smear, MultipleScattering);
    if (simTrack.GetHits().size() < (unsigned int) MinimumNumberOfHits) {
      LOG(6, "Hit has only " << simTrack.GetHits().size() << " hits, skip it");
      continue;
    }
    // Add the true track to the visualization
    gVis.AddTrack(simTrack);
    // only check measurement model if no detector resolution simulation and
    // no MS effects
    //  check_measurement_model(t, simTrack.GetHits()); 
    LOG(5, "Track " << i << " has " << simTrack.GetHits().size() << " hits.");
    // show hits only for the first 100 tracks
    if (i < 100) {
      for (vector<Hit *>::const_iterator it = simTrack.GetHits().begin(); 
	   it != simTrack.GetHits().end(); it++) {
	gVis.AddHit(*(*it));
      }
    }
    // now that we have a simulated track and hits, let us reconstruct a track
    // from the measured hits.
    // cout << "True track parameters: " << simTrack.GetGlobalParameters();
    try {
      /// reconstruct track using nominal coordinate system
      LOG(5, "reco track, ");
      Track recoTrack = tracking.MakeRecoTrack(simTrack.GetHits());
//       vis.AddTrack(recoTrack);
      // // consistency check: in case without misalignment and ideal
      // // hits (i.e. no detector resolution), the reconstructed track should
      // // match the true track ideally.
      // if ((recoTrack.GetGlobalParameters()-simTrack.GetGlobalParameters()).norm() > 1E-8) {
      // 	cout << endl << "bad tracking" << endl;
      // 	cout << "true track" << simTrack.GetGlobalParameters() << endl
      // 	     << "refence track" << recoTrack.GetGlobalParameters() << endl;
      // }

      // fill histograms only if more information is requested
      if (gLogLevel > 2) {
	fill_tracking_pulls(recoTrack, simTrack);
	fill_track_chi2(recoTrack);
      }
 
      //////////////////////////////////////////////////////////////////////
      // fill event information
      myEvent->SetEvent(i);
      myEvent->SetTrackNumber(1);
      KalmanFilterAlignmentInputProvider kfaip(recoTrack);
      kfaip.Fill(*myEvent);
      // fill tree
      AlignTree->Fill();
    }
    CATCH;
    // cout << "True track parameters: " << simTrack.GetGlobalParameters();
    i++;
  }

  // Show detectors, tracks & hits
  gVis.Show();

  gBenchmark->Show("simulation");

  // Data I/O
  alignment_data->Write();
  alignment_data->Close();
  delete alignment_data;

  // Run interactive interface
  theApp->Run();
  delete theApp;

  delete gDetector;

  return 0;
}
