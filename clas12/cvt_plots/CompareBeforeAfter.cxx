#include <cstdlib>
#include <list>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TStyle.h>
#include <TLine.h>
#include <TText.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TArc.h>
#include "../../event/AlignEvent.h"


#include <iostream>



#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"
#include "TGraphErrors.h"


using namespace std;

TString formatRvsR(char* x, char* y){
  return Form("%s vs %s;Residual %s [mm];Residual %s [mm]", x,y, x,y);
  
}

//change mode to "count" to make it show the number of hits on a given module
TH2F* drawPositions(int i1, int i2){
  gPad->SetMargin(.2,.15,.15, .1);
  double R = 235;
  int Nbins = 235;
  TH2F*  hxy =   new TH2F(Form("temp %d %d",i1,i2), ";x (mm);y (mm);", Nbins,-R,R,Nbins,-R,R);
  hxy->Draw("COLZ");
  //cout << prof->GetNbinsX() << "bins" << endl;
  for(int mm = 0; mm<102; mm++){
    if(mm < 84){
      int m = mm%42;
      
      double pi = TMath::Pi();
      double phim =  m<10 ? 2*pi*m/10 : (m<24 ? 2*pi*(m-10)/14 : 2*pi*(m-24)/18);
      double dphim = m<10 ? pi/10 : (m<24 ? pi/14 : pi/18);
      //double rm = m<10 ? 65 : (m<24 ? 93 : 120) + 2.7*(mm>=42);
      double rm = (m<10 ? 65 : (m<24 ? 93 : 120));
      rm += 4*(mm>=42);
      phim -= pi/2;
      while(phim < -pi)
        phim += 2*pi;
      while(phim > pi)
        phim -= 2*pi;
      TLine* line = new TLine();
      if (mm == i1)
        line->SetLineColor(kRed);
      else if (mm == i2)
        line->SetLineColor(kBlue);
      else
        line->SetLineColor(kGray);
      line->DrawLine(rm*cos(phim+dphim)/cos(dphim), rm*sin(phim+dphim)/cos(dphim),rm*cos(phim-dphim)/cos(dphim), rm*sin(phim-dphim)/cos(dphim));
    }
      /*for(int i = 0; i<hxy->GetNbinsX(); i++){
        double x = -R +2*i*R/Nbins;
        
        for(int j = 0; j<hxy->GetNbinsY(); j++){
          double y = -R +2*j*R/Nbins;
          double phi = TMath::ATan2(y,x);
          double rc = x*cos(phim)+y*sin(phim);
          if(mm<42 && rc>rm && rc<rm+2.3 && abs(phi-phim)<dphim){
            hxy->SetBinContent(i,j,n);
          }
          else if(mm>=42 && rc>rm+2.7 && rc<rm+5.0 && abs(phi-phim)<dphim){
            hxy->SetBinContent(i,j,n);
          }
       }
      }*/
    else if(mm>=84){ //BMT:
      int m = mm-84;

      int layer = m/3;
      int sector = m%3;

      double radii[] = {146,161,176,191,206,221};
      double rm = radii[layer], drm = 4;
      double phisector = -TMath::Pi()/6+sector*TMath::Pi()*2/3;
      double dphi = TMath::Pi()/3-0.1;//the width between sectors is not right, but whatever.
      TArc* arc = new TArc(0,0, radii[layer], (phisector-dphi)*180/TMath::Pi(),(phisector+dphi)*180/TMath::Pi());
      if (mm == i1)
        arc->SetLineColor(kRed);
      else if (mm == i2)
        arc->SetLineColor(kBlue);
      else
        arc->SetLineColor(kGray);
      arc->SetFillStyle(0);
      arc->SetNoEdges(1);
      arc->Draw();
      //arc->DrawArc(0,0, radii[layer], (phisector-dphi)*180/TMath::Pi(),(phisector+dphi)*180/TMath::Pi());
      //cout<<  mm << " " << layer << " " << sector << " " << n << " " << phisector*180/TMath::Pi() << endl;
    }
  }
  hxy->SetStats(0);
  
  return hxy;
}

double dmin=-30,dmax = 30;
double phimin=-TMath::Pi(), phimax = TMath::Pi();
double tandipmin= -1,tandipmax =1.5;
double zmin=-70, zmax =10;
TString plotsDir ="";

TF1* g = new TF1("gaus","gaus", -.1,.1);
double getSigma(TH1* h){
  double mu = h->GetMean();
  double sigma = h->GetStdDev();
  
  double minfit = mu-1.5*sigma;
  double maxfit = mu+1.5*sigma;
  g->SetParameter("Mean", mu);
  g->SetParameter("Sigma", sigma);
  h->Fit(g, "N", "", minfit, maxfit);
  mu = g->GetParameter("Mean");
  sigma = g->GetParameter("Sigma");
  return sigma;
}

//calculates a quantile-truncated standard deviation
double getStdTruncated(TH1* h, double q=0.10){
  TH1* cum = h->GetCumulative();
  cum->Scale(1/h->Integral());
  TH1* clo = (TH1*)h->Clone();
  for(int i = 0; i< cum->GetNbinsX(); i++){
    if (cum->GetBinContent(i)<q/2 || cum->GetBinContent(i)>1-q/2)
      {
	clo->SetBinContent(i,0);
      }
  }
  return clo->GetStdDev();
}

double getMeanTruncated(TH1* h, double q=0.10){
  if (q==0){
    return h->GetMean();
  }
   TH1* cum = h->GetCumulative();
  cum->Scale(1/h->Integral());
  TH1* clo = (TH1*)h->Clone();
  for(int i = 0; i< cum->GetNbinsX(); i++){
    if (cum->GetBinContent(i)<q/2 || cum->GetBinContent(i)>1-q/2)
      {
        clo->SetBinContent(i,0);
      }
  }
  return clo->GetMean();
}

double getMode(TH1* h){
  int maxbin = h->GetMaximumBin();
  double max = h->GetBinContent(maxbin);
  double left = h->GetXaxis()->GetXmin();
  double right = h->GetXaxis()->GetXmax();
  double prev = 0;
  for(int j = maxbin;j>0; j--){
    double bc = h->GetBinContent(j);
    if(bc<=max/2 && prev>= max/2 && j<maxbin){
      left = (h->GetXaxis()->GetBinCenter(j)+h->GetXaxis()->GetBinCenter(j-1))/2;
      break;
    }
    prev=bc;
  }
  prev = 0;
  for (int j = maxbin; j<h->GetNbinsX(); j++){
    double bc = h->GetBinContent(j);
    if(bc<=max/2 && prev>= max/2 && j>maxbin){
      right = (h->GetXaxis()->GetBinCenter(j)+h->GetXaxis()->GetBinCenter(j-1))/2;
      break;
    }
    prev = bc;
  }
  double val= h->GetBinCenter(maxbin);
  if (abs(val) <0.0001)
    val=0;
  return val;
}


double getFWHM(TH1* h){
  int maxbin = h->GetMaximumBin();
  double max = h->GetBinContent(maxbin);
  double left = h->GetXaxis()->GetXmin();
  double right = h->GetXaxis()->GetXmax();
  double prev = 0;
  for(int j = maxbin;j>0; j--){
    double bc = h->GetBinContent(j);
    if(bc<=max/2 && prev>= max/2 && j<maxbin){
      left = (h->GetXaxis()->GetBinCenter(j)+h->GetXaxis()->GetBinCenter(j-1))/2;
      break;
    }
    prev=bc;
  }
  prev = 0;
  for (int j = maxbin; j<h->GetNbinsX(); j++){
    double bc = h->GetBinContent(j);
    if(bc<=max/2 && prev>= max/2 && j>maxbin){
      right = (h->GetXaxis()->GetBinCenter(j)+h->GetXaxis()->GetBinCenter(j-1))/2;
      break;
    }
    prev = bc;
  }
  return right-left;
}




TGraphErrors* createProfile(TH2D * h, int color, int markerstyle, double shift,double ywindow=1.5, int fitType=1){

  int n =h->GetXaxis()->GetNbins();
  double x[n];
  double y[n];
  double ex[n];
  double ey[n];

  double largestMu = 0;
  double largestMuError = 0;
  double chi2 = 0;
  for(int i = 0; i<n;i++){
    TH1* proj = h->ProjectionY("temp", i,i+1);
    
    if (fitType==2){ //full width at half max
      int maxbin = proj->GetMaximumBin();                                                                                                             
      double max = proj->GetBinContent(maxbin);
      double left = proj->GetXaxis()->GetXmin();
      double right = proj->GetXaxis()->GetXmax();
      double prev = 0;
      for(int j = maxbin;j>0; j--){
        double bc = proj->GetBinContent(j);
        if(bc<=max/2 && prev>= max/2 && j<maxbin){
          left = (proj->GetXaxis()->GetBinCenter(j)+proj->GetXaxis()->GetBinCenter(j-1))/2;
	  break;
	}
	prev=bc;
      }
      prev = 0;
      for (int j = maxbin; j<proj->GetNbinsX(); j++){
	double bc = proj->GetBinContent(j);
        if(bc<=max/2 && prev>= max/2 && j>maxbin){
          right = (proj->GetXaxis()->GetBinCenter(j)+proj->GetXaxis()->GetBinCenter(j-1))/2;
          break;
        }
	prev = bc;
      }
      x[i] =  h->GetXaxis()->GetBinCenter(i+1)+shift;
      //x[i] = h->GetMean();
      ex[i] = 0;

      double yi=proj->GetBinCenter(maxbin);
      
      /*if(proj->GetBinContent(maxbin+1)>proj->GetBinContent(maxbin-1)){
	yi=(proj->GetXaxis()->GetBinCenter(maxbin)*proj->GetBinContent(maxbin)+proj->GetXaxis()->GetBinCenter(maxbin+1)*proj->GetBinContent(maxbin+1))/
	  (proj->GetBinContent(maxbin)+proj->GetBinContent(maxbin+1));
      }
      if(proj->GetBinContent(maxbin+1)<proj->GetBinContent(maxbin-1)){
        yi=(proj->GetXaxis()->GetBinCenter(maxbin)*proj->GetBinContent(maxbin)+proj->GetXaxis()->GetBinCenter(maxbin-1)*proj->GetBinContent(maxbin-1))/
          (proj->GetBinContent(maxbin)+proj->GetBinContent(maxbin-1));
      }*/ 
       
      y[i] = proj->GetMean();
      cout << x[i] << " " << y[i] << " " << yi << endl;
      ey[i] = (right-left)/2;
      cout << "FWHM" << x[i] << " " << y[i] << " " << ey[i] << endl;
      if(proj->GetEntries()<10){
        y[i]=100;
      }
      continue;
    }
    
    double mu = proj->GetMean();
    //double dmu = proj->GetMeanError();
    double sigma = proj->GetStdDev();
    
    double minfit = mu-2*sigma;
    double maxfit = mu+2*sigma;
    g->SetParameter("Mean", mu);
    g->SetParameter("Sigma", sigma);
    proj->Fit(g, "NQ", "", minfit, maxfit);
    mu = g->GetParameter("Mean");
    sigma = g->GetParameter("Sigma");
    double dmu = g->GetParError(1);
    if(sigma>5 || fitType == 0){
      mu = getMeanTruncated(proj,0);
      sigma = getStdTruncated(proj,0);
      //proj->Draw();
      //gPad->SaveAs(plotsDir+ "/error.pdf");
      
    }
    x[i] = h->GetXaxis()->GetBinCenter(i+1)+shift;
    ex[i] = 0;
    y[i] = mu;
    ey[i] = sigma;
    if(abs(mu)>largestMu){
      largestMu = abs(mu);
      largestMuError = dmu;
      
    }
    //for empty or nearly empty bins, hide the marker
    if(proj->GetEntries()<10){
      y[i]=100;
    }
    chi2+=mu*mu/(dmu*dmu);
  }
  gStyle->SetTitleFontSize(0.08);
  TGraphErrors * profile = new TGraphErrors(n, x, y, ex, ey);
  profile->SetMarkerStyle(markerstyle);
//profile->SetMarkerSize(1);
  profile->SetMarkerColor(color);
  profile->SetLineColor(color);
  profile->SetLineStyle(0);
  profile->GetXaxis()->SetRangeUser(h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  profile->SetTitle(h->GetTitle());
  profile->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  profile->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
  profile->GetYaxis()->SetTitleOffset(h->GetYaxis()->GetTitleOffset());
  profile->GetHistogram()->SetTitleSize(h->GetTitleSize(), "XYZT");
  profile->GetHistogram()->SetLabelSize(h->GetLabelSize(), "XYZT");
  profile->GetXaxis()->SetTitleSize(h->GetTitleSize("X"));
  profile->GetYaxis()->SetTitleSize(h->GetTitleSize("Y"));
  profile->GetXaxis()->SetLabelSize(h->GetLabelSize("X"));
  profile->GetYaxis()->SetLabelSize(h->GetLabelSize("Y"));
  
  profile->GetHistogram()->SetMaximum(ywindow);   // along
  profile->GetHistogram()->SetMinimum(-ywindow);  //   Y
  
  //cout << h->GetName() << " worst mu: " << largestMu << "+-"<< largestMuError << endl;
  //cout << h->GetName() << " chi2: " << chi2 << endl;
  return profile;
//return h->ProfileX(h->GetName()+(TString&)"_prof");
  
}

int main(int argc, char * argv[]) {
  // Record start time
  //auto start = std::chrono::high_resolution_clock::now();
  
  TString inputFileNameBefore;
  TString inputFileNameAfter;
  //TString plotsDir;

  double maxResid=5;
  vector<TString> inputFiles;
  TString label ="";
  bool isMC = 0;
  int fitType = 1;
  double scaleWindows=1;
  for(Int_t i=1;i<argc;i++){
    TString opt=argv[i];
    if((opt.Contains("--before="))){
      inputFileNameBefore=opt(9,opt.Sizeof());
    }
    if((opt.Contains("--scaleWindows="))){
      scaleWindows=((TString)opt(15,opt.Sizeof())).Atof();
    }
    if((opt.Contains("--after="))){
      inputFileNameAfter=opt(8,opt.Sizeof());
    }
    if (opt.Contains("--plotsdir=")){
      plotsDir = opt(11,opt.Sizeof());
    }
    if (opt.Contains("--maxResid=")){
      maxResid=((TString)opt(11,opt.Sizeof())).Atof();
    }
    if (opt.Contains("-l")){
      i++;
      label = argv[i];
    }
    if (opt == "--isMC"){
      isMC=1;
    }
    if (opt == "--meanStd"){
      fitType = 0;
    }
    if (opt == "--fwhm"){
      fitType = 2;
      //exit(0);
    }
  }
  
  

  gStyle->SetTitleSize(0.06,"XYZT");
  gStyle->SetLabelSize(0.06,"XYZT");
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadBottomMargin(.13);
  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);

  //gStyle->SetTitleSize(0.05,"XYZT");
  //gStyle->SetLabelSize(0.05,"XYZT");
  gStyle->SetPadLeftMargin(.10);
  gStyle->SetPadRightMargin(.075);
  TCanvas* c2 = new TCanvas("c2","c2",1200,900);

  gStyle->SetTitleSize(0.06,"XYZT");
  gStyle->SetLabelSize(0.06,"XYZT");
  gStyle->SetPadLeftMargin(.19);
  gStyle->SetPadRightMargin(.03);
  gStyle->SetPadBottomMargin(.13);
  //gStyle->SetPadTopMargin(0.1);
  TCanvas* c3 = new TCanvas("c3","c3",1200,1600);
  c3->Divide(3,4);
  
  TLegend* legend1 = new TLegend(0.153, 0.7, 0.977, 0.898);
  TLegend* legend2 = new TLegend(0.153, 0.7, 0.977, 0.898);
  TLegend* legend3 = new TLegend(0.153, 0.7, 0.977, 0.898);
  TLegend* legend4 = new TLegend(0.153, 0.7, 0.70, 0.898);
  legend1->SetBorderSize(0);
  legend2->SetBorderSize(0);
  legend3->SetBorderSize(0);
  legend4->SetBorderSize(0);
  
  TLegend* legend5 = new TLegend(0.13, 0.74, 0.40, 0.89);
  TLegend* legend6 = new TLegend(0.25, 0.67, 0.7, 0.87);
  legend5->SetBorderSize(0);
  legend6->SetBorderSize(0);

  
  gStyle->SetPadBottomMargin(.12);
  gStyle->SetPadLeftMargin(0.13);
  TCanvas* c4 = new TCanvas("c4","c4",800,1300);
  c4->Divide(3,5);
  
  TCanvas* c5 = new TCanvas("c5","c5",800,1300);
  c5->Divide(3,5);
  
  
  for(int before_after = 0; before_after<2; before_after++){
    
    int n_chi2_lt_2 = 0;
    int n_tracks = 0;
    // auto save every MB
    //AlignTree->SetAutoSave(1000000);
    
    // these histograms are for debugging.  The limits and labels for
    // some of these histograms are taylored for the Clas12 CVT.
    
    /*if(before_after){
      gStyle->SetMarkerColor(kBlack);
      gStyle->SetLineColor(kBlack);
    } else {
      gStyle->SetMarkerColor(kRed);
      gStyle->SetLineColor(kRed);
    }*/
    
    TString suffix = before_after ? "After" : "Before";
    gStyle->SetTitleSize(0.05, "XYZT");
    gStyle->SetLabelSize(0.05, "XYZT");
    TH1F*  hchi2ndof = new TH1F("hchi2ndof"+suffix, "track #chi^{2}/ndof;track #chi^{2}/n_{dof};# of tracks", 100, 0, 20);
    
    TH1F* residuals_svt = new TH1F ("res_svt"+suffix, "SVT residuals;residual [mm];# of clusters", 500, -1.0*scaleWindows, 1.0*scaleWindows);
    TH1F* residuals_bmtz = new TH1F ("res_bmtz"+suffix, "BMTZ residuals;residual [mm];# of clusters", 500, -4.0*scaleWindows, 4.0*scaleWindows);
    TH1F* residuals_bmtc = new TH1F ("res_bmtc"+suffix, "BMTC residuals;residual [mm];# of clusters", 500, -2.0*scaleWindows, 2.0*scaleWindows);

    cout << "500 bins" << endl;
    
    hchi2ndof->GetYaxis()->SetMaxDigits(3);
    residuals_svt->GetYaxis()->SetMaxDigits(3);
    residuals_bmtz->GetYaxis()->SetMaxDigits(3);
    residuals_bmtc->GetYaxis()->SetMaxDigits(3);
    
    //shift the "before" by a fraction of a bin for clarity
    double shift_module = -(1.-before_after)*1./5;


    

    TString statlbl = " (gaus fit #mu#pm#sigma)";
    if(fitType==2) statlbl = " (mean#pm FWHM/2)";
    if(fitType==0) statlbl = " (mean#pm std)";

    gStyle->SetTitleSize(0.06, "XYZT");
    gStyle->SetLabelSize(0.06, "XYZT");
    #define RESID_BINS 101, -3,3
    TH2D* residuals_vs_module =  new TH2D ("res_mod_"+suffix, "Residuals (all modules);module # ;residual [mm]"+statlbl, 102, 0+shift_module, 102+shift_module,
                                           RESID_BINS);
    
    //sets the error bar to be the std dev, instead of standard error on the mean
    //residuals_vs_module->SetErrorOption("s");
    
    int color = before_after ? kBlack: kRed;
     residuals_svt->SetLineColor(color);
     residuals_bmtz->SetLineColor(color);
     residuals_bmtc->SetLineColor(color);
     hchi2ndof->SetLineColor(color);
    
    residuals_vs_module->SetLineColor(color);
    residuals_vs_module->SetMarkerColor(color);
    
     int linestyle = before_after ? 1: 2;
     residuals_svt->SetLineStyle(linestyle);
     residuals_bmtz->SetLineStyle(linestyle);
     residuals_bmtc->SetLineStyle(linestyle);
     hchi2ndof->SetLineStyle(linestyle);
    
    int markerstyle=before_after ? 21: 25;
    residuals_vs_module->SetMarkerStyle(markerstyle);
    
    
    //number of bins for the residuals vs kinematics plots.
    int nbins_kin = 15;
    //shift the "before" markers by a 5th of a bin for visual clarity
    double shift_d  = -(1.-before_after)*(dmax-dmin)/nbins_kin/5;
    double shift_phi  = -(1.-before_after)*(phimax-phimin)/nbins_kin/5;
    double shift_z  = -(1.-before_after)*(zmax-zmin)/nbins_kin/5;
    double shift_tandip  = -(1.-before_after)*(tandipmax-tandipmin)/nbins_kin/5;
    cout << shift_d << endl;
    
    
    gStyle->SetLabelSize(0.06,"XYZT");
    gStyle->SetTitleSize(0.06, "XYZT");
    TH2D* residuals_vs_phi_svt = new  TH2D ("res_phi_svt"+suffix, "SVT residuals vs #phi_{0};#phi_{0} [rad];residual [mm]"+ statlbl, nbins_kin,phimin, phimax, RESID_BINS); residuals_vs_phi_svt->SetTitleSize(0.1, "XYZT"); residuals_vs_phi_svt->SetLabelSize(0.06, "XYZT");
    TH2D* residuals_vs_d0_svt =  new TH2D ("res_d0_svt"+suffix, "SVT residuals vs d_{0};d_{0} [mm]; residual [mm]" + statlbl, nbins_kin, dmin, dmax, RESID_BINS); residuals_vs_d0_svt->SetTitleSize(0.06, "XYZT"); residuals_vs_d0_svt->SetLabelSize(0.06, "XYZT");
    TH2D* residuals_vs_theta_svt =  new TH2D ("res_theta_svt"+suffix, "SVT residuals vs t_{0};t_{0};residual [mm]"+statlbl, nbins_kin, tandipmin, tandipmax, RESID_BINS); residuals_vs_theta_svt->SetTitleSize(0.06, "XYZT"); residuals_vs_theta_svt->SetLabelSize(0.06, "XYZT");
    TH2D* residuals_vs_z_svt =  new TH2D ("res_z_svt"+suffix, "SVT residuals vs z_{0};z_{0} [mm];residual [mm]"+statlbl, nbins_kin, zmin, zmax, RESID_BINS);  residuals_vs_z_svt->SetTitleSize(0.06, "XYZT"); residuals_vs_z_svt->SetLabelSize(0.06, "XYZT");
    
    TH2D* residuals_vs_phi_bmtz = new  TH2D ("res_phi_bmtz"+suffix, "BMTZ residuals vs #phi_{0};#phi_{0} [rad];residual [mm]"+statlbl, nbins_kin, phimin,phimax, RESID_BINS); residuals_vs_phi_bmtz->SetTitleSize(0.06, "XYZT"); residuals_vs_phi_bmtz->SetLabelSize(0.06, "XYZT");
    TH2D* residuals_vs_d0_bmtz = new  TH2D ("res_d0_bmtz"+suffix, "BMTZ residuals vs d_{0};d_{0} [mm]; residual [mm]"+statlbl, nbins_kin, dmin, dmax, RESID_BINS); residuals_vs_d0_bmtz->SetTitleSize(0.06, "XYZT"); residuals_vs_d0_bmtz->SetLabelSize(0.06, "XYZT\
");
    TH2D* residuals_vs_theta_bmtz = new  TH2D ("res_theta_bmtz"+suffix, "BMTZ residuals vs t_{0};t_{0};residual [mm]"+statlbl, nbins_kin, tandipmin, tandipmax, RESID_BINS); residuals_vs_theta_bmtz->SetTitleSize(0.06, "XYZT"); residuals_vs_theta_bmtz->SetLabelSize(0.06, "XYZT\
");
    TH2D* residuals_vs_z_bmtz =  new TH2D ("res_z_bmtz"+suffix, "BMTZ residuals vs z_{0};z_{0} [mm];residual [mm]"+statlbl, nbins_kin, zmin, zmax, RESID_BINS); residuals_vs_z_bmtz->SetTitleSize(0.06, "XYZT"); residuals_vs_z_bmtz->SetLabelSize(0.06, "XYZT\
");
    
    TH2D* residuals_vs_phi_bmtc = new  TH2D ("res_phi_bmtc"+suffix, "BMTC residuals vs #phi_{0};#phi_{0} [rad];residual [mm]"+statlbl, nbins_kin, phimin, phimax,RESID_BINS);
    TH2D* residuals_vs_d0_bmtc =  new TH2D ("res_d0_bmtc"+suffix, "BMTC residuals vs d_{0};d_{0} [mm]; residual [mm]"+statlbl, nbins_kin, dmin, dmax, RESID_BINS);
    TH2D* residuals_vs_theta_bmtc =  new TH2D ("res_theta_bmtc"+suffix, "BMTC residuals vs t_{0};t_{0};residual [mm]"+statlbl, nbins_kin, tandipmin, tandipmax, RESID_BINS);
    TH2D* residuals_vs_z_bmtc =  new TH2D ("res_z_bmtc"+suffix, "BMTC residuals vs z_{0};z_{0} [mm];residual [mm]"+statlbl, nbins_kin, zmin, zmax, RESID_BINS);
    
/*#define format_profile(profile,min,max) \
    
    
    format_profile(residuals_vs_d0_svt,-0.6,0.6)
    format_profile(residuals_vs_phi_svt,-0.6,0.6)
    format_profile(residuals_vs_z_svt,-0.6,0.6)
    format_profile(residuals_vs_theta_svt,-0.6,.6)
    format_profile(residuals_vs_d0_bmtz,-5,5)
    format_profile(residuals_vs_phi_bmtz,-5,5)
    format_profile(residuals_vs_z_bmtz,-5,5)
    format_profile(residuals_vs_theta_bmtz,-5,5)
    format_profile(residuals_vs_d0_bmtc,-5,5)
    format_profile(residuals_vs_phi_bmtc,-5,5)
    format_profile(residuals_vs_z_bmtc,-5,5)
    format_profile(residuals_vs_theta_bmtc,-5,5)
*/
    //int fontSizePrev = gStyle->GetTitleFontSize();
    //gStyle->SetTitleFontSize(0.12);
    gStyle->SetTitleSize(0.06,"XYZT");
    gStyle->SetLabelSize(0.06,"XYZT");
    
    vector<int> module1_resid_vs_resid;
    vector<int> module2_resid_vs_resid;
    vector<TH2 *> h_module_resid_vs_resid;
    int nbinsRvR = 50;
    double maxResidRvR = 1.0;
    //double maxResidRvR = before_after? 0.6 : 2.4;
    int nDivisions = 604;
    //back to back SVT
    module1_resid_vs_resid.push_back(5);
    module2_resid_vs_resid.push_back(47);
    h_module_resid_vs_resid.push_back(new TH2F("hrvsr1"+suffix, suffix+": "+formatRvsR("SVT L1S6","SVT L2S6"), nbinsRvR, -maxResidRvR, maxResidRvR, nbinsRvR, -maxResidRvR, maxResidRvR));
    
    //opposite sides of the detector (cosmics only)
    module1_resid_vs_resid.push_back(6);
    module2_resid_vs_resid.push_back(1);
    h_module_resid_vs_resid.push_back(new TH2F("hrvsr2"+suffix, formatRvsR("SVT L1S7","SVT L1S2"), nbinsRvR, -maxResidRvR, maxResidRvR, nbinsRvR, -maxResidRvR, maxResidRvR));
    //overlapping SVT sectors
    module1_resid_vs_resid.push_back(5);
    module2_resid_vs_resid.push_back(9+10+14);
    h_module_resid_vs_resid.push_back(new TH2F("hrvsr3"+suffix, formatRvsR("SVT L1S6","SVT L5S10"), nbinsRvR, -maxResidRvR, maxResidRvR, nbinsRvR, -maxResidRvR, maxResidRvR));
    //svt vs bmtc
    module1_resid_vs_resid.push_back(5);
    module2_resid_vs_resid.push_back(85);
    h_module_resid_vs_resid.push_back(new TH2F("hrvsr4"+suffix, formatRvsR("SVT L1S6","BMTC L1S2"), nbinsRvR, -maxResidRvR, maxResidRvR, nbinsRvR, -maxResidRvR, maxResidRvR));
    //svt vs bmtz
    maxResidRvR=3;
    module1_resid_vs_resid.push_back(5);
    module2_resid_vs_resid.push_back(88);
    h_module_resid_vs_resid.push_back(new TH2F("hrvsr5"+suffix, formatRvsR("SVT L1S6","BMTZ L2S2"), nbinsRvR, -maxResidRvR, maxResidRvR, nbinsRvR, -maxResidRvR, maxResidRvR));
    //bmtz vs bmtz (same sector)
    module1_resid_vs_resid.push_back(88);
    module2_resid_vs_resid.push_back(97);
    h_module_resid_vs_resid.push_back(new TH2F("hrvsr6"+suffix, suffix+": "+formatRvsR("BMTZ L2S2","BMTZ L5S2"), nbinsRvR, -maxResidRvR, maxResidRvR, nbinsRvR, -maxResidRvR, maxResidRvR));
    //bmtz vs bmtz (different sector)
    module1_resid_vs_resid.push_back(88);
    module2_resid_vs_resid.push_back(96);
    h_module_resid_vs_resid.push_back(new TH2F("hrvsr7"+suffix, formatRvsR("BMTZ L2S2","BMTZ L5S1"), nbinsRvR, -maxResidRvR, maxResidRvR, nbinsRvR, -maxResidRvR, maxResidRvR));
    //bmtz vs bmtc (same sector)
    module1_resid_vs_resid.push_back(88);
    module2_resid_vs_resid.push_back(100);
    h_module_resid_vs_resid.push_back(new TH2F("hrvsr8"+suffix, formatRvsR("BMTZ L2S2","BMTC L6S2"), nbinsRvR, -maxResidRvR, maxResidRvR, nbinsRvR, -maxResidRvR, maxResidRvR));
    //bmtc vs bmtc (same sector)
    maxResidRvR = 1;
    module1_resid_vs_resid.push_back(85);
    module2_resid_vs_resid.push_back(100);
    h_module_resid_vs_resid.push_back(new TH2F("hrvsr9"+suffix, formatRvsR("BMTC L1S2","BMTC L6S2"), nbinsRvR, -maxResidRvR, maxResidRvR, nbinsRvR, -maxResidRvR, maxResidRvR));
    //bmtc vs bmtc (diff. sector)
    module1_resid_vs_resid.push_back(85);
    module2_resid_vs_resid.push_back(99);
    h_module_resid_vs_resid.push_back(new TH2F("hrvsr10", formatRvsR("BMTZ L1S2","BMTC L6S1"), nbinsRvR, -maxResidRvR, maxResidRvR, nbinsRvR, -maxResidRvR, maxResidRvR));
    cout << "check" << endl;
    for (int i = 0; i<h_module_resid_vs_resid.size();i++){
      TH2* h = h_module_resid_vs_resid[i];
      h->SetTitleSize(0.06,"XYZT");
      h->SetLabelSize(0.06,"XYZT");
      h->GetXaxis()->SetNdivisions(nDivisions);
      h->GetYaxis()->SetNdivisions(nDivisions);
      //h->SetPadLeftMargin(0.14);
      //h->SetPadBottomMargin(0.14);
    }
    
    for (int i =0; i<10; i++){
      if(i<5)
        c4->cd(3*i+1);
      else
        c5->cd(3*(i-5)+1);
      cout <<module1_resid_vs_resid[i] << "  " << module2_resid_vs_resid[i] << endl;
      TH1*h = drawPositions(module1_resid_vs_resid[i],module2_resid_vs_resid[i]);
      h->SetTitleSize(0.06,"XYZT");
      h->SetLabelSize(0.06,"XYZ");
      h->GetXaxis()->SetNdivisions(405);
      h->GetXaxis()->SetNdivisions(405);
      //h->Draw("COLZ");
      
    }
    
    
    //gStyle->SetTitleFontSize(fontSizePrev);
    
    cout << "initialized histograms" << endl;
    int events = 0, tracks=0;
    
    int ndof =0;
    float chi2 =0;
    
    //if there is no input file
    
    TString inputFileName = before_after ? inputFileNameAfter : inputFileNameBefore;
    if(inputFileName==TString())  {
      std::cout << " *** please provide a file name..." << std::endl;
      exit(0);
    }
    
    
    
    /////////////////////////////////////
    
    TFile * inputFile = TFile::Open(inputFileName);
    cout << "found file"<<endl;
    
    AlignInfo * alignInfo = (AlignInfo *)inputFile->Get("AlignInfo");
    if(alignInfo == NULL){
      cout << "align info not found" << endl;
      exit(0);
    }
    cout << "align info"  << endl;
    int alignables = alignInfo->GetNAlignables();
    int parameters = alignInfo->GetNParameters();
    cout << alignables << "alignables; " <<  parameters << " parameters"<< endl;
    
    AlignEvent *aevent = new AlignEvent();
    TTree *AlignTree = (TTree*)inputFile->Get("AlignTree");
    AlignTree->SetBranchAddress("AlignEvent", & aevent);
    
    float d0;
    float phi;
    float z;
    float tandip;
    AlignTree->SetBranchAddress("q0", &d0);
    AlignTree->SetBranchAddress("q1", &phi);
    AlignTree->SetBranchAddress("q2", &z);
    AlignTree->SetBranchAddress("q3", &tandip);
    
    
    for(int i = 0; i<AlignTree->GetEntriesFast();i++){
      
      AlignTree->GetEntry(i);
      int ndof = aevent->GetNdof();
      float chi2 = aevent->GetChi2();
      
      //cout << d0 << " "<< phi << " "<< z << " " << tandip << endl;
      
      hchi2ndof->Fill(chi2/ndof);
      if(chi2/ndof<2)
        n_chi2_lt_2++;
      n_tracks++;
      
      
      for(int j = 0; j < aevent->GetMeasuredCovariance()->GetNrows(); j++){
        double resid = (*aevent->GetMeasurements())(j)-(*aevent->GetTrackPrediction())(j);
        double res = sqrt((*aevent->GetMeasuredCovariance())[j][j]);
        double module = aevent->GetIndex()->At(j);
        
        
        residuals_vs_module->Fill(module,resid);
        
        //residual vs residual plots
        for (int jj =0;jj<module1_resid_vs_resid.size(); jj++){
          if (module != module1_resid_vs_resid[jj])
            continue;
          for (int k =0;k<aevent->GetMeasuredCovariance()->GetNrows(); k++){
              int module2 = aevent->GetIndex()->At(k);
            if (module2 != module2_resid_vs_resid[jj])
              continue;
            double resid2 = (*aevent->GetMeasurements())(k)-(*aevent->GetTrackPrediction())(k);
            h_module_resid_vs_resid[jj]->Fill(resid,resid2);
          }
        }
        
       
        if(module >=84 && module < 102){
          if (module <=86 || module>=93 && module<=95 || module>=99){
            residuals_bmtc->Fill(resid);
            residuals_vs_d0_bmtc->Fill(d0, resid);
            residuals_vs_phi_bmtc->Fill(phi, resid);
            residuals_vs_z_bmtc->Fill(z, resid);
            residuals_vs_theta_bmtc->Fill(tandip, resid);
          } else {
            residuals_bmtz->Fill(resid);
            residuals_vs_d0_bmtz->Fill(d0, resid);
            residuals_vs_phi_bmtz->Fill(phi, resid);
            residuals_vs_z_bmtz->Fill(z, resid);
            residuals_vs_theta_bmtz->Fill(tandip, resid);
          }
        } else if (module == 102){
          //beamspot.
        } else {
          residuals_svt->Fill(resid);
          residuals_vs_d0_svt->Fill(d0, resid);
          residuals_vs_phi_svt->Fill(phi, resid);
          residuals_vs_z_svt->Fill(z, resid);
          residuals_vs_theta_svt->Fill(tandip, resid);
        }
        
      }
      
    }
    cout << suffix << endl;
    
    cout << suffix << ":  frac of events with chi2/ndof<2: "  << (n_chi2_lt_2/(double)n_tracks) << endl;
    
    TString opt = before_after ? "SAME" : "";
    
    if (isMC){
      residuals_svt->SetTitle((TString)"    " + residuals_svt->GetTitle()+(TString)" (MC)");
      residuals_bmtz->SetTitle((TString)"    " + residuals_bmtz->GetTitle()+(TString)" (MC)");
      residuals_bmtc->SetTitle((TString)"    " + residuals_bmtc->GetTitle()+(TString)" (MC)");
      hchi2ndof->SetTitle((TString)"    " + hchi2ndof->GetTitle()+(TString)" (MC)");
    }
    
    c1->cd(1);
    //legend1->AddEntry(residuals_svt, Form(suffix + ",\n RMS = %.2f mm, fit #sigma = %.3f mm",residuals_svt->GetRMS(), getSigma(residuals_svt)),"l");
    legend1->AddEntry(residuals_svt, Form(suffix + ",\n mean = %.0f #mum, FWHM = %.0f #mum",1000*getMeanTruncated(residuals_svt,0), 1000*getFWHM(residuals_svt),"l"));
      //legend1->AddEntry(residuals_svt, Form(suffix + ",\n mode = %.0f #mum, FWHM = %.0f #mum",1000*getMode(residuals_svt), 1000*getFWHM(residuals_svt),"l")); 
    residuals_svt->SetTitleSize(0.05, "T");  
    residuals_svt->Draw(opt);
    residuals_svt->SetMaximum(residuals_svt->GetMaximum()*(isMC ? 7 : 5.5));
    c1->cd(2);
    //legend2->AddEntry(residuals_bmtz, Form(suffix + ",\n RMS = %.2f mm, fit #sigma = %.2f mm",residuals_bmtz->GetRMS(), getSigma(residuals_bmtz)),"l");
    legend2->AddEntry(residuals_bmtz, Form(suffix + ",\n mean = %.0f #mum, FWHM = %.0f #mum",1000*getMeanTruncated(residuals_bmtz,0), 1000*getFWHM(residuals_bmtz),"l"));
    residuals_bmtz->SetTitleSize(0.05, "T");
    residuals_bmtz->Draw(opt);
    residuals_bmtz->SetMaximum(residuals_bmtz->GetMaximum()*(isMC ? 9 : 7));
    c1->cd(3);
    //legend3->AddEntry(residuals_bmtc, Form(suffix+ ",\n RMS = %.2f mm, fit #sigma = %.2f mm",residuals_bmtc->GetRMS(), getSigma(residuals_bmtc)), "l");
    legend3->AddEntry(residuals_bmtc, Form(suffix+ ",\n mean = %.0f #mum, FWHM = %.0f #mum",1000*getMeanTruncated(residuals_bmtc,0), 1000*getFWHM(residuals_bmtc), "l"));
    residuals_bmtc->SetTitleSize(0.05, "T");
    residuals_bmtc->Draw(opt);
    residuals_bmtc->SetMaximum(residuals_bmtc->GetMaximum()*(isMC ? 3 : 3.5));
    c1->cd(4);
    legend4->AddEntry(hchi2ndof, Form(suffix+", mean = %.1f",hchi2ndof->GetMean()),"l");
    hchi2ndof->SetMaximum(hchi2ndof->GetMaximum()*(isMC ? 21 : 16));
    hchi2ndof->SetTitleSize(0.05, "T");
    hchi2ndof->Draw(opt);

   
    
    opt = before_after ? "SP" : "AP";
    
    TLine* line = new TLine();
    line->SetLineStyle(2);
    line->SetLineColorAlpha(kBlack,0.5);
    c2->cd();
    residuals_vs_module->GetYaxis()->SetTitleOffset(1.5);
    residuals_vs_module->SetLabelSize(0.05, "XYZT");
    residuals_vs_module->SetTitleSize(0.05, "XYZT");
    line->DrawLine(0,0,102,0);
    //residuals_vs_module->SetMinimum(-2);
    //residuals_vs_module->SetMaximum(3);
    if (isMC){
      suffix = "MC (" + suffix + ")";
      residuals_vs_module->SetTitle("Residuals (all modules, MC)");
    }
    TGraphErrors* t = createProfile(residuals_vs_module, color, markerstyle, shift_module,1.5,fitType);
    if(before_after){
      double worstMuSVT=0;
      double worstMuBMTZ=0;
      double worstMuBMTC=0;
      for (int k = 0; k<t->GetN(); k++){
        double mu = t->GetPointY(k);
        if(k >=84 && k < 102){
          if (k <=86 || k>=93 && k<=95 || k>=99){
            if (abs(mu)>worstMuBMTC)
              worstMuBMTC=abs(mu);
          } else {
            if (abs(mu)>worstMuBMTZ){
              worstMuBMTZ=abs(mu);
            }
          }
        } else {
          if(abs(mu)>worstMuSVT){
            worstMuSVT=abs(mu);
          }
        }
      }
      cout << "worst residuals for modules (SVT, BMTZ, BMTC): " << worstMuSVT <<", " << worstMuBMTZ << ", "<< worstMuBMTC<<  endl;
    }
    t->GetHistogram()->SetMinimum(-1.15);
    t->GetHistogram()->SetMaximum(1.7);
    t->GetHistogram()->GetYaxis()->SetTitleOffset(1);
    t->Draw(opt);
    
    legend5->AddEntry(residuals_vs_module,suffix, "lp");
    
    if (isMC){
      residuals_vs_d0_svt->SetTitle(residuals_vs_d0_svt->GetTitle()+(TString)" (MC)");
      residuals_vs_d0_bmtz->SetTitle(residuals_vs_d0_bmtz->GetTitle()+(TString)" (MC)");
      residuals_vs_d0_bmtc->SetTitle(residuals_vs_d0_bmtc->GetTitle()+(TString)" (MC)");
      
      residuals_vs_phi_svt->SetTitle(residuals_vs_phi_svt->GetTitle()+(TString)" (MC)");
      residuals_vs_phi_bmtz->SetTitle(residuals_vs_phi_bmtz->GetTitle()+(TString)" (MC)");
      residuals_vs_phi_bmtc->SetTitle(residuals_vs_phi_bmtc->GetTitle()+(TString)" (MC)");
      
      residuals_vs_z_svt->SetTitle(residuals_vs_z_svt->GetTitle()+(TString)" (MC)");
      residuals_vs_z_bmtz->SetTitle(residuals_vs_z_bmtz->GetTitle()+(TString)" (MC)");
      residuals_vs_z_bmtc->SetTitle(residuals_vs_z_bmtc->GetTitle()+(TString)" (MC)");
      
      residuals_vs_theta_svt->SetTitle(residuals_vs_theta_svt->GetTitle()+(TString)" (MC)");
      residuals_vs_theta_bmtz->SetTitle(residuals_vs_theta_bmtz->GetTitle()+(TString)" (MC)");
      residuals_vs_theta_bmtc->SetTitle(residuals_vs_theta_bmtc->GetTitle()+(TString)" (MC)");
    }
    
    //gStyle->SetTitleH(0.06);
    c3->cd(1);TGraphErrors *graph = createProfile(residuals_vs_d0_svt, color, markerstyle, shift_d,0.5,fitType); graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06,"XYZT"); graph->Draw(opt);
    legend6->AddEntry(graph,suffix, "lp"); graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06,"XYZT");
    line->DrawLine(dmin,0,dmax,0);
    c3->cd(4); graph=createProfile(residuals_vs_phi_svt, color, markerstyle, shift_phi,0.5,fitType);graph->Draw(opt); graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06, "XYZT");
    line->DrawLine(phimin,0,phimax,0);
    c3->cd(7); graph=createProfile(residuals_vs_z_svt, color, markerstyle, shift_z,0.5,fitType);graph->Draw(opt); graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06, "XYZT");
    line->DrawLine(zmin,0,zmax,0);
    c3->cd(10); graph=createProfile(residuals_vs_theta_svt, color, markerstyle, shift_tandip,0.5,fitType);graph->Draw(opt); graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06, "XYZT");
    line->DrawLine(tandipmin,0,tandipmax,0);
    
    
    c3->cd(2); graph=createProfile(residuals_vs_d0_bmtz, color, markerstyle, shift_d, 1.5,fitType);graph->Draw(opt);graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06, "XYZT");
    line->DrawLine(dmin,0,dmax,0);
    c3->cd(5); graph=createProfile(residuals_vs_phi_bmtz, color, markerstyle, shift_phi, 1.5,fitType);graph->Draw(opt); graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06, "XYZT");
    line->DrawLine(phimin,0,phimax,0);
    c3->cd(8); graph=createProfile(residuals_vs_z_bmtz, color, markerstyle, shift_z, 1.5,fitType);graph->Draw(opt); graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06, "XYZT");
    line->DrawLine(zmin,0,zmax,0);
    c3->cd(11); graph=createProfile(residuals_vs_theta_bmtz, color, markerstyle, shift_tandip, 1.5,fitType);graph->Draw(opt); graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06, "XYZT");
    line->DrawLine(tandipmin,0,tandipmax,0);
    
    c3->cd(3); graph=createProfile(residuals_vs_d0_bmtc, color, markerstyle, shift_d, 1.5,fitType);graph->Draw(opt); graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06, "XYZT");
    line->DrawLine(dmin,0,dmax,0);
    c3->cd(6); graph=createProfile(residuals_vs_phi_bmtc, color, markerstyle, shift_phi, 1.5,fitType);graph->Draw(opt); graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06, "XYZT");
    line->DrawLine(phimin,0,phimax,0);
    c3->cd(9); graph=createProfile(residuals_vs_z_bmtc, color, markerstyle, shift_z, 1.5,fitType);graph->Draw(opt); graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06, "XYZT");
    line->DrawLine(zmin,0,zmax,0);
    c3->cd(12); graph=createProfile(residuals_vs_theta_bmtc, color, markerstyle, shift_tandip, 1.5,fitType);graph->Draw(opt); graph->GetHistogram()->SetLabelSize(0.06, "XYZT"); graph->GetHistogram()->SetTitleSize(0.06, "XYZT");
    //c3->cd(12);residuals_vs_theta_bmtc->Draw(opt);
    line->DrawLine(tandipmin,0,tandipmax,0);
    
    for (int i = 0; i<10; i++){
      if (i<5){
        c4->cd(2+3*i+before_after);
      }
      else{
        c5->cd(2+3*(i-5)+before_after);
      }
      gPad->SetMargin(.2,.15,.15, .1);
     
      gStyle->SetOptStat(0);
      h_module_resid_vs_resid[i]->Draw("COLZ");
      //TPaveText *pt = (TPaveText*)(gPad->GetPrimitive("title"));
      //cout << pt<<endl;
      //cout << (h_module_resid_vs_resid[i]->GetListOfPrimitives()->At(0)->GetName()) << endl;
      //pt->SetTextSize(0.07);
    }
    
    
    inputFile->Close();
  }
  
   TString tag =  "";
    if (fitType== 0)
      tag =  "_mean_std";
    if (fitType ==2)
      tag = "_FWHM";
  TString ext = "pdf";
  if(!(plotsDir==TString())){
    c1->cd(1);
    legend1->Draw();
    c1->cd(2);
    legend2->Draw();
    c1->cd(3);
    legend3->Draw();
    c1->cd(4);
    legend4->Draw();
    c1->SaveAs(plotsDir+"/residuals_1d"+tag+"." +ext);
    
    c2->cd();
    legend5->Draw();
    TLine* line = new TLine();
    line->SetLineStyle(3);
    double xs[] = {42,84,87,93,96,99};
    for(double x : xs){
      if(x == 84)
        line->SetLineWidth(2);
      else
        line->SetLineWidth(1);
      line->DrawLine(x,-1.15,x,x <= 84? 1.7: 1.4);
    }
    
    TText *text = new TText();
    text->SetTextSize(.06);
    text->DrawText(12, 0.82, "SVT (inner)");
    text->DrawText(12+42, 0.82, "SVT (outer)");
    
    text->DrawText(90, 1.5, "BMT");
    text->SetTextSize(.04);
    
    text->DrawText(84+0.8, 1.3, "C");
    text->DrawText(87+0.8, 1.3, "Z");
    text->DrawText(90+0.8, 1.3, "Z");
    text->DrawText(93+0.8, 1.3, "C");
    text->DrawText(96+0.8, 1.3, "Z");
    text->DrawText(99+0.8, 1.3, "C");
  
    cout << "tag is " << tag << ", fitType=" << fitType <<  endl;
    c2->SaveAs(plotsDir+"/residuals_module"+tag+ "." + ext);
    c3->cd(1);
    legend6->Draw();
    //if(label != ""):
    //  text->DrawText(0, label);
    c3->SaveAs(plotsDir+"/residuals_kinematics"+ tag+"." + ext);
    gStyle->SetTitleFontSize(0.07);
    
    c4->SaveAs(plotsDir+"/rvr1."+ext);
    c5->SaveAs(plotsDir+"/rvr2."+ext);
  }
  
  
  
  //auto finish = std::chrono::high_resolution_clock::now();
  // std::chrono::duration<double> elapsed = finish - start;
  //std::cout << "Elapsed time: " << elapsed.count()<< "s, events = "<<events<< "\n";
  
}
